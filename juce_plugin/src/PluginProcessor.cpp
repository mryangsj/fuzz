#include "PluginProcessor.h"
#include "PluginEditor.h"

//==============================================================================
AudioPluginAudioProcessor::AudioPluginAudioProcessor()
     : AudioProcessor (BusesProperties()
                     #if ! JucePlugin_IsMidiEffect
                      #if ! JucePlugin_IsSynth
                       .withInput  ("Input",  juce::AudioChannelSet::stereo(), true)
                      #endif
                       .withOutput ("Output", juce::AudioChannelSet::stereo(), true)
                     #endif
                       ),
       SZ_h5File ("/Users/yangshijie/Desktop/fuzz/SZ.hdf5", HighFive::File::ReadOnly)
{
    //==============================================================================
    // Grab the parameter with the identifier from the APVTS and cast it to an
    // AudioParameter object, so that we can quickly lookup the value in the latter.
    castParameter(apvts, MyParameterID::input_gain, inputGainParameter);
    castParameter(apvts, MyParameterID::lin, linParameter);
    castParameter(apvts, MyParameterID::vol, volParameter);

    //==============================================================================
    // Listening to parameter changes
    apvts.state.addListener(this);

    //==============================================================================
    // Load the neural network model
    loadNN();


    //==============================================================================
}

AudioPluginAudioProcessor::~AudioPluginAudioProcessor()
{
    apvts.state.removeListener(this);
}

//==============================================================================
const juce::String AudioPluginAudioProcessor::getName() const
{
    return JucePlugin_Name;
}

bool AudioPluginAudioProcessor::acceptsMidi() const
{
   #if JucePlugin_WantsMidiInput
    return true;
   #else
    return false;
   #endif
}

bool AudioPluginAudioProcessor::producesMidi() const
{
   #if JucePlugin_ProducesMidiOutput
    return true;
   #else
    return false;
   #endif
}

bool AudioPluginAudioProcessor::isMidiEffect() const
{
   #if JucePlugin_IsMidiEffect
    return true;
   #else
    return false;
   #endif
}

double AudioPluginAudioProcessor::getTailLengthSeconds() const
{
    return 0.0;
}

int AudioPluginAudioProcessor::getNumPrograms()
{
    return 1;   // NB: some hosts don't cope very well if you tell them there are 0 programs,
                // so this should be at least 1, even if you're not really implementing programs.
}

int AudioPluginAudioProcessor::getCurrentProgram()
{
    return 0;
}

void AudioPluginAudioProcessor::setCurrentProgram (int index)
{
    juce::ignoreUnused (index);
}

const juce::String AudioPluginAudioProcessor::getProgramName (int index)
{
    juce::ignoreUnused (index);
    return {};
}

void AudioPluginAudioProcessor::changeProgramName (int index, const juce::String& newName)
{
    juce::ignoreUnused (index, newName);
}

//==============================================================================
void AudioPluginAudioProcessor::prepareToPlay (double sampleRate, int samplesPerBlock)
{
    // Use this method as the place to do any pre-playback
    // initialisation that you need..
    juce::ignoreUnused (samplesPerBlock);

    //==============================================================================
    // Update the sample rate and related parameters
    if (float(Fs) != float(sampleRate)) { set_fs(sampleRate); }

    // This forces updateParameters() to be executed the very first time processBlock
    // is called.
    refreshParameters();
    parametersChanged.store(true);

    //==============================================================================
    // Set initial conditions
    setInitalConditions();
}

void AudioPluginAudioProcessor::releaseResources()
{
    // When playback stops, you can use this as an opportunity to free up any
    // spare memory, etc.
}

bool AudioPluginAudioProcessor::isBusesLayoutSupported (const BusesLayout& layouts) const
{
  #if JucePlugin_IsMidiEffect
    juce::ignoreUnused (layouts);
    return true;
  #else
    // This is the place where you check if the layout is supported.
    // In this template code we only support mono or stereo.
    // Some plugin hosts, such as certain GarageBand versions, will only
    // load plugins that support stereo bus layouts.
    if (layouts.getMainOutputChannelSet() != juce::AudioChannelSet::mono()
     && layouts.getMainOutputChannelSet() != juce::AudioChannelSet::stereo())
        return false;

    // This checks if the input layout matches the output layout
   #if ! JucePlugin_IsSynth
    if (layouts.getMainOutputChannelSet() != layouts.getMainInputChannelSet())
        return false;
   #endif

    return true;
  #endif
}

void AudioPluginAudioProcessor::processBlock (juce::AudioBuffer<float>& buffer,
                                              juce::MidiBuffer& midiMessages)
{
    juce::ignoreUnused (midiMessages);

    juce::ScopedNoDenormals noDenormals;
    auto totalNumInputChannels  = getTotalNumInputChannels();
    auto totalNumOutputChannels = getTotalNumOutputChannels();

    // In case we have more outputs than inputs, this code clears any output
    // channels that didn't contain input data, (because these aren't
    // guaranteed to be empty - they may contain garbage).
    // This is here to avoid people getting screaming feedback
    // when they first compile a plugin, but obviously you don't need to keep
    // this code if your algorithm always overwrites all the output channels.
    for (auto i = totalNumInputChannels; i < totalNumOutputChannels; ++i)
        buffer.clear (i, 0, buffer.getNumSamples());


    //==============================================================================
    // These two lines of code do a thread-safe check to see whether
    // "parametersChanged" is true. If so, it calls "updateParameters()" method to
    // perform the parameter calculations. It also immediately sets
    // "parametersChanged" back to false.
    // This "check if true and set back to false" is a single atomic operation.
    // It either succeeds or it fails. If the operation fails, some other thread was
    // trying to write to "parametersChanged" while "processBlock()" was trying to
    // read it. This is not a big deal. You'll miss the parameter update in this
    // block, but most likely the operation will succeed in one of the next blocks
    // and "updateParameters()" will eventually get called.
    bool expected = true;
    if (parametersChanged.compare_exchange_strong(expected, false)) { updateParameters(); }


    // This is the place where you'd normally do the guts of your plugin's
    // audio processing...
    // Make sure to reset the state if your inner loop is processing
    // the samples and the outer loop is handling the channels.
    // Alternatively, you can process the samples with the channels
    // interleaved by keeping the same state.
    for (int channel = 0; channel < totalNumInputChannels; ++channel) {
        auto* channelData = buffer.getWritePointer (channel);
        auto ch = (size_t)channel;

        for (int i = 0; i < buffer.getNumSamples(); ++i) {
            //==============================================================================
            // Manage memory elements states
            C71_b[ch] = C71_a[ch];
            C81_b[ch] = C81_a[ch];
            C91_b[ch] = C91_a[ch];
            C94_b[ch] = C94_a[ch];


            //==============================================================================
            // Forward scanning
            S71_a(1) = C71_b[ch];
            S71_a(2) = inputAmp * channelData[i];
            S71_b(0) = S71_S.row(0) * S71_a;

            P72_a(2) = S71_b(0);
            P72_b(0) = P72_S.row(0) * P72_a;

            P81_a(1) = C81_b[ch];
            P81_b(0) = P81_S.row(0) * P81_a;

            S91_a(2) = C91_b[ch];
            S91_b(0) = S91_S.row(0) * S91_a;

            P92_a(2) = S91_b(0);
            P92_b(0) = P92_S.row(0) * P92_a;

            S93_a(2) = P92_b(0);
            S93_b(0) = S93_S.row(0) * S93_a;

            S94_a(1) = C94_b[ch];
            S94_a(2) = S93_b(0);
            S94_b(0) = S94_S.row(0) * S94_a;

            //==============================================================================
            // R-type junction scattering
            RJ_af << P72_b(0),
                     P81_b(0),
                     S94_b(0),
                     Vcc;

            RJ_bf = RJ_Sf[lin_index][vol_index] * RJ_af;


            //==============================================================================
            // Local scattering
            BJT_in << RJ_bf,
                      RJ_Zr[lin_index][vol_index];

            // Neural network forward scanning
            gn1.forward(BJT_in);
            linear1.forward(gn1.outs);
            relu1.forward(linear1.outs);
            linear2.forward(relu1.outs);
            relu2.forward(linear2.outs);
            linear3.forward(relu2.outs + relu1.outs);


            //==============================================================================
            // Backward scanning
            RJ_ab << linear3.outs,
                     RJ_af;

            RJ_bb = RJ_Sb[lin_index][vol_index] * RJ_ab;

            P72_a(0) = RJ_bb(0);
            P72_b(2) = P72_S.row(2) * P72_a;

            S71_a(0) = P72_b(2);
            S71_b(1) = S71_S.row(1) * S71_a;

            P81_a(0) = RJ_bb(1);
            P81_b(1) = P81_S.row(1) * P81_a;

            S94_a(0) = RJ_bb(2);
            S94_b(1) = S94_S.row(1)* S94_a;
            S94_b(2) = S94_S.row(2) * S94_a;

            S93_a(0) = S94_b(2);
            S93_b(2) = S93_S.row(2) * S93_a;

            P92_a(0) = S93_b(2);
            P92_b(2) = P92_S.row(2) * P92_a;

            S91_a(0) = P92_b(2);
            S91_b(1) = S91_S.row(1) * S91_a;
            S91_b(2) = S91_S.row(2) * S91_a;


            //==============================================================================
            // Update memory elements states
            C71_a[ch] = S71_b(1);
            C81_a[ch] = P81_b(1);
            C91_a[ch] = S91_b(2);
            C94_a[ch] = S94_b(1);


            //==============================================================================
            // Collect initial conditions
//            float eps = 1e-7f;
//            if (
//                    std::abs(C71_b[ch] - C71_a[ch]) < eps
//                    && std::abs(C81_b[ch] - C81_a[ch]) < eps
//                    && std::abs(C91_b[ch] - C91_a[ch]) < eps
//                    && std::abs(C94_b[ch] - C94_a[ch]) < eps
//                    ) {
//                std::cout << "C71_ic = \n" << std::fixed << std::setprecision(32) << C71_a[ch] << std::endl;
//                std::cout << "C81_ic = \n" << std::fixed << std::setprecision(32) << C81_a[ch] << std::endl;
//                std::cout << "C91_ic = \n" << std::fixed << std::setprecision(32) << C91_a[ch] << std::endl;
//                std::cout << "C94_ic = \n" << std::fixed << std::setprecision(32) << C94_a[ch] << std::endl;
//            }


            //==============================================================================
            // Output
            channelData[i] = (S91_b(1) / 2.0f) / 4.0f;
//            DBG((S91_b(1) / 2.0f) / 4.0f);
        }
    }
}

//==============================================================================
bool AudioPluginAudioProcessor::hasEditor() const
{
    return true; // (change this to false if you choose to not supply an editor)
}

juce::AudioProcessorEditor* AudioPluginAudioProcessor::createEditor()
{
//    return new AudioPluginAudioProcessorEditor (*this);

    auto editor = new juce::GenericAudioProcessorEditor(*this);
//    editor->setSize(400, 300);
    return editor;
}

//==============================================================================
void AudioPluginAudioProcessor::getStateInformation (juce::MemoryBlock& destData)
{
    // You should use this method to store your parameters in the memory block.
    // You could do that either as raw data, or use the XML or ValueTree classes
    // as intermediaries to make it easy to save and load complex data.
//    juce::ignoreUnused (destData);

    // Save plugin's current state and load it back the next time the plug-in is used.
    copyXmlToBinary(*apvts.copyState().createXml(), destData);
    //DBG(apvts.copyState().toXmlString());
}

void AudioPluginAudioProcessor::setStateInformation (const void* data, int sizeInBytes)
{
    // You should use this method to restore your parameters from this memory block,
    // whose contents will have been created by the getStateInformation() call.
//    juce::ignoreUnused (data, sizeInBytes);


    // You are provided with a block of binary data. First, this calls
    // "getXmlFromBinary()" to parse the binary data into an XML document. Then it
    // verifies this XML contains the <Parameters> tag. And finally, it calls
    // "apvts.replaceState()" to update the values of all the parameters to the new
    // values.

    // After loading and restoring the plugin's state, you set the
    // "parametersChanged" boolean to true to signal "processBlock()" that it should
    // call "updateParameters()" again.

    // There is no way to know exactly when "setStateInformation()" may be called in
    // the lifetime of the plug-in, so you should assume it could happen at any time,
    // even if the plug-in is busy processing audio. Hence, it's necessary to
    // recalculate anything that depends on these parameter values the very next
    // time "processBlock()" is called.

    // The APVTS "replaceState()" method is thread-safe, by the way. That's good
    // thing too, because there is no guarantee that "getStateInformation()" or
    // "setStateInformation()" will be called from any particular thread, it doesn't
    // have to be the UI thread.
    std::unique_ptr<juce::XmlElement> xml(getXmlFromBinary(data, sizeInBytes));
    if (xml != nullptr && xml->hasTagName(apvts.state.getType())) {
        apvts.replaceState(juce::ValueTree::fromXml(*xml));
        parametersChanged.store(true);
    }
}

//==============================================================================
// This creates new instances of the plugin..
juce::AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
    return new AudioPluginAudioProcessor();
}

//==============================================================================
juce::AudioProcessorValueTreeState::ParameterLayout AudioPluginAudioProcessor::createParameterLayout()
{
    //==============================================================================
    juce::AudioProcessorValueTreeState::ParameterLayout layout;


    //==============================================================================
    // add "INPUT_GAIN" parameter
    const juce::NormalisableRange<float> range_input_gain(-30.0f, 30.0f, 1e-1f);
    layout.add(std::make_unique<juce::AudioParameterFloat>(
            MyParameterID::input_gain,
            "INPUT GAIN",
            range_input_gain,
            0.0f,
            juce::AudioParameterFloatAttributes().withLabel("dB")));


    //==============================================================================
    // add "LIN" parameter
    const juce::NormalisableRange<float> range_lin(0.0f, 100.0f, 1e-0f);
    layout.add(std::make_unique<juce::AudioParameterFloat>(
            MyParameterID::lin,
            "LIN",
            range_lin,
            50.0f,
            juce::AudioParameterFloatAttributes().withLabel("%")));


    //==============================================================================
    // add "VOL" parameter
    const juce::NormalisableRange<float> range_vol(0.0f, 100.0f, 1e-0f);
    layout.add(std::make_unique<juce::AudioParameterFloat>(
            MyParameterID::vol,
            "LEVEL",
            range_vol,
            50.0f,
            juce::AudioParameterFloatAttributes().withLabel("%")));

    //==============================================================================
    return layout;
}

//==============================================================================
void AudioPluginAudioProcessor::refreshParameters()
{
    inputGain = inputGainParameter->get();
    inputAmp = juce::Decibels::decibelsToGain(inputGain);

    alpha_lin = linParameter->get() / 100.0f;
    lin_index = (size_t)(std::round(alpha_lin * 100.0f));

    alpha_vol = volParameter->get() / 100.0f;
    vol_index = (size_t)(std::round(alpha_vol * 100.0f));

//    pDiode = diodeParameter->getIndex();
//    pConfig = configParameter->getIndex();
//    diode = static_cast<DiodeType>(pDiode);
//    config = static_cast<DiodeConfig>(pConfig);
//    diodeClipper.set_model(diode, config);
}

void AudioPluginAudioProcessor::updateParameters()
{
    if (lastChangedParamID == MyParameterID::input_gain.getParamID()) {
        inputGain = inputGainParameter->get();
        inputAmp = juce::Decibels::decibelsToGain(inputGain);
    } else if (lastChangedParamID == MyParameterID::lin.getParamID()) {
        alpha_lin = linParameter->get() / 100.0f;
        lin_index = (size_t)(std::round(alpha_lin * 100.0f));
    } else if (lastChangedParamID == MyParameterID::vol.getParamID()) {
        alpha_vol = volParameter->get() / 100.0f;
        vol_index = (size_t)(std::round(alpha_vol * 100.0f));
        if (vol_index == 0) { alpha_vol = 1e-4f; }
        update_port9_scattering_matrix(Fs, alpha_vol);
    }
}

void AudioPluginAudioProcessor::setInitalConditions() {
    for (size_t i = 0; i < 2; ++i) {
        C71_a[i] = C71_ic;
        C81_a[i] = C81_ic;
        C91_a[i] = C91_ic;
        C94_a[i] = C94_ic;
    }
}

//==============================================================================
void AudioPluginAudioProcessor::loadNN() {
    std::string hdf5Path = "/Users/yangshijie/Desktop/fuzz/models/2x16_gn.h5";
    HighFive::File hdf5File = HighFive::File(hdf5Path, HighFive::File::ReadOnly);

    loadBatchNorm1d(hdf5File, "gn1", gn1);
    loadDense(hdf5File, "linear1", linear1);
    loadDense(hdf5File, "linear2", linear2);
    loadDense(hdf5File, "linear3", linear3);
}

//==============================================================================
template <typename BatchNormType>
void AudioPluginAudioProcessor::loadBatchNorm1d(const HighFive::File& h5File, const std::string& layerName, BatchNormType& batchNorm1d) {
    auto runningMeanDataset = h5File.getDataSet(layerName + "/running_mean");
    auto runningVarianceDataset = h5File.getDataSet(layerName + "/running_var");
    auto epsilonDataset = h5File.getDataSet(layerName + "/eps");
    auto weightsDataset = h5File.getDataSet(layerName + "/weight");
    auto biasDataset = h5File.getDataSet(layerName + "/bias");

    auto runningMean = runningMeanDataset.read<std::vector<float>>();
    auto runningVariance = runningVarianceDataset.read<std::vector<float>>();
    auto epsilon = epsilonDataset.read<float>();
    auto weights = weightsDataset.read<std::vector<float>>();
    auto bias = biasDataset.read<std::vector<float>>();

    runningMeanDataset.read(runningMean);
    runningVarianceDataset.read(runningVariance);
    epsilonDataset.read(epsilon);
    weightsDataset.read(weights);
    biasDataset.read(bias);

    batchNorm1d.setRunningMean(runningMean);
    batchNorm1d.setRunningVariance(runningVariance);
    batchNorm1d.setEpsilon(epsilon);
    batchNorm1d.setGamma(weights);
    batchNorm1d.setBeta(bias);
}

//==============================================================================
template <typename DenseType>
void AudioPluginAudioProcessor::loadDense(const HighFive::File& h5File, const std::string& layerName, DenseType& dense)
{
    auto weightsDataset = h5File.getDataSet(layerName + "/weight");
    auto biasDataset = h5File.getDataSet(layerName + "/bias");

    auto weights = weightsDataset.read<std::vector<std::vector<float>>>();
    auto bias = biasDataset.read<std::vector<float>>();

    weightsDataset.read(weights);
    biasDataset.read(bias);

    dense.setWeights(weights);
    dense.setBias(bias.data());
}

//==============================================================================
void AudioPluginAudioProcessor::set_fs(double sampleRate) {
    Fs = float(sampleRate);
    T = float(1 / sampleRate);

    update_port7_scattering_matrix(Fs);
    update_port8_scattering_matrix(Fs);
    update_port9_scattering_matrix(Fs, alpha_vol);
    load_RJ_scattering_matrix(Fs);
}

//==============================================================================
void AudioPluginAudioProcessor::
update_port7_scattering_matrix(float fs) {
    // Update impedance of fs-dependent elements
    C71_z = 1.0f / (2.0f * fs * C71);

    // Update impedance vectors
    S71_Zv << C71_z + R71_z,
              C71_z,
              R71_z;

    P72_Zv << S71_Zv(0) * R72_z / (S71_Zv(0) + R72_z),
              R72_z,
              S71_Zv(0);

    // Recompute scattering matrices
    S71_S = get_scattering_matrix_series<3>(S71_Bv, S71_Zv);
    P72_S = get_scattering_matrix_parallel<3>(P72_Qv, P72_Zv);
}

void AudioPluginAudioProcessor::
update_port8_scattering_matrix(float fs) {
    // Update impedance of fs-dependent elements
    C81_z = 1.0f / (2.0f * fs * C81);

    // Update impedance vectors
    P81_Zv << C81_z * R81_z / (C81_z + R81_z),
              C81_z,
              R81_z;

    // Recompute scattering matrices
    P81_S = get_scattering_matrix_parallel<3>(P81_Qv, P81_Zv);
}

void AudioPluginAudioProcessor::
update_port9_scattering_matrix(float fs, float alpha_vol_) {
    // Update impedance of user parameters dependent elements
    C91_z = 1.0f / (2.0f * fs * C91);
    C94_z = 1.0f / (2.0f * fs * C94);

    R92_z = alpha_vol_ * Rvol;
    R93_z = (1 - alpha_vol_) * Rvol;

    // Update impedance vectors
    S91_Zv << R91_z + C91_z,
              R91_z,
              C91_z;

    P92_Zv << R92_z * S91_Zv(0) / (R92_z + S91_Zv(0)),
              R92_z,
              S91_Zv(0);

    S93_Zv << R93_z + P92_Zv(0),
              R93_z,
              P92_Zv(0);

    S94_Zv << C94_z + S93_Zv(0),
              C94_z,
              S93_Zv(0);

    // Recompute scattering matrices
    S91_S = get_scattering_matrix_series<3>(S91_Bv, S91_Zv);
    P92_S = get_scattering_matrix_parallel<3>(P92_Qv, P92_Zv);
    S93_S = get_scattering_matrix_series<3>(S93_Bv, S93_Zv);
    S94_S = get_scattering_matrix_series<3>(S94_Bv, S94_Zv);
}

//==============================================================================
template<size_t NumPortsTotal>
Eigen::Matrix<float, NumPortsTotal, NumPortsTotal> AudioPluginAudioProcessor::
get_scattering_matrix_series(
        const Eigen::Vector<float, NumPortsTotal>& Bv,
        const Eigen::Vector<float, NumPortsTotal>& Zv) const {
    const auto B = Bv.transpose();
    const auto Z = Zv.asDiagonal();
    const auto I = Eigen::Matrix<float, NumPortsTotal, NumPortsTotal>::Identity();
    const auto ZBt = Z * B.transpose();
    return I - 2.0f * ZBt * (B * ZBt).inverse() * B;
}

template<size_t NumPortsTotal>
Eigen::Matrix<float, NumPortsTotal, NumPortsTotal> AudioPluginAudioProcessor::
get_scattering_matrix_parallel(
        const Eigen::Vector<float, NumPortsTotal>& Qv,
        const Eigen::Vector<float, NumPortsTotal>& Zv) const {
    const auto Q = Qv.transpose();
    const auto Z = Zv.asDiagonal();
    const auto I = Eigen::Matrix<float, NumPortsTotal, NumPortsTotal>::Identity();
    const auto QZ_inv = Q * Z.inverse();
    return 2.0f * Q.transpose() * (QZ_inv * Q.transpose()).inverse() * QZ_inv - I;
}

//==============================================================================
bool AudioPluginAudioProcessor::
load_RJ_scattering_matrix(float fs) {
    std::string fs_str {};
    try {
        switch (fs_dict.at(fs)) {
            case FS::RATE_44100:
                fs_str = "44.1k";
                break;
            case FS::RATE_48000:
                fs_str = "48.0k";
                break;
            case FS::RATE_88200:
                fs_str = "88.2k";
                break;
            case FS::RATE_96000:
                fs_str = "96.0k";
                break;
            case FS::RATE_176400:
                fs_str = "176.4k";
                break;
            case FS::RATE_192000:
                fs_str = "192.0k";
                break;
        }
    } catch (const std::out_of_range& e) {
        DBG("Unsupported sampling rate detected: " + std::to_string(fs));
        return false;
    }

    for (size_t i = 0; i < 101; ++i) {
        for (size_t j = 0; j < 101; ++j) {
            std::string Sf_path = "/Sf/" + fs_str + "/lin=" + std::to_string(i) + "%_vol=" + std::to_string(j) + "%";
            std::string Sb_path = "/Sb/" + fs_str + "/lin=" + std::to_string(i) + "%_vol=" + std::to_string(j) + "%";
            std::string Zr_path = "/Zr/" + fs_str + "/lin=" + std::to_string(i) + "%_vol=" + std::to_string(j) + "%";

            auto Sf_vector = H5Easy::load<std::vector<std::vector<float>>>(SZ_h5File, Sf_path);
            auto Sb_vector = H5Easy::load<std::vector<std::vector<float>>>(SZ_h5File, Sb_path);
            auto Zr_vector = H5Easy::load<std::vector<std::vector<float>>>(SZ_h5File, Zr_path);

            for (size_t ii = 0; ii < 2; ++ii) {
                for (size_t jj = 0; jj < 4; ++jj) {
                    RJ_Sf[i][j](ii, jj) = Sf_vector[ii][jj];
                }
            }

            for (size_t ii = 0; ii < 3; ++ii) {
                for (size_t jj = 0; jj < 6; ++jj) {
                    RJ_Sb[i][j](ii, jj) = Sb_vector[ii][jj];
                }
            }

            RJ_Zr[i][j](0) = Zr_vector[0][0];
            RJ_Zr[i][j](1) = Zr_vector[0][1];
            RJ_Zr[i][j](2) = Zr_vector[1][1];
        }
    }

    return true;
}
#pragma once

#include <juce_audio_processors/juce_audio_processors.h>
#include <RTNeural/RTNeural.h>
#include <highfive/H5Easy.hpp>
#include <highfive/highfive.hpp>
//#include <highfive/eigen.hpp>
#include <chrono>


namespace MyParameterID
{
    #define PARAMETER_ID(str) const juce::ParameterID str(#str, 1);

    PARAMETER_ID(input_gain)
    PARAMETER_ID(lin)
    PARAMETER_ID(vol)

    #undef PARAMETER_ID
}


//==============================================================================
class AudioPluginAudioProcessor final : public juce::AudioProcessor,
                                        private juce::ValueTree::Listener
{
public:
    //==============================================================================
    AudioPluginAudioProcessor();
    ~AudioPluginAudioProcessor() override;

    //==============================================================================
    void prepareToPlay (double sampleRate, int samplesPerBlock) override;
    void releaseResources() override;
    bool isBusesLayoutSupported (const BusesLayout& layouts) const override;
    void processBlock (juce::AudioBuffer<float>&, juce::MidiBuffer&) override;
    using AudioProcessor::processBlock;

    //==============================================================================
    juce::AudioProcessorEditor* createEditor() override;
    bool hasEditor() const override;

    //==============================================================================
    const juce::String getName() const override;

    bool acceptsMidi() const override;
    bool producesMidi() const override;
    bool isMidiEffect() const override;
    double getTailLengthSeconds() const override;

    //==============================================================================
    int getNumPrograms() override;
    int getCurrentProgram() override;
    void setCurrentProgram (int index) override;
    const juce::String getProgramName (int index) override;
    void changeProgramName (int index, const juce::String& newName) override;

    //==============================================================================
    void getStateInformation (juce::MemoryBlock& destData) override;
    void setStateInformation (const void* data, int sizeInBytes) override;

private:
    //==============================================================================
    // JUCE has a powerful class called "ValueTree", a tree structure that can
    // contain objects of any type. The APVTS uses such a "ValueTree" to hold the
    // plug-in's parameters.
    // One of its useful features is that the "ValueTree" allows you to connect
    // listeners to nodes in the tree, for receiving notifications when something in
    // the tree changes. It can also serialize the tree to XML, supports undo/redo
    // functionality, automatically manages the lifetimes of nodes through reference
    // counting and more. It's usually added to the public section, so that the
    // editor class that handles the UI can also access it.
    juce::AudioProcessorValueTreeState apvts {
        *this,
        nullptr,
        "PARAMETERS",
        createParameterLayout()
    };


    //==============================================================================
    // "createParameterLayout()" will be called by APVTS. Inside this method is
    // where you will instantiate all the AudioParameter objects.
    static juce::AudioProcessorValueTreeState::ParameterLayout createParameterLayout();


    //==============================================================================
    // This method is used to grab the parameter with the identifier from APVTS, and
    // casts it to the destination AudioParameter object.
    template<typename T>
    inline static void castParameter(
            juce::AudioProcessorValueTreeState& apvts,
            const juce::ParameterID& id,
            T& destination) {
        destination = dynamic_cast<T>(apvts.getParameter(id.getParamID()));
        jassert(destination); // parameter does not exist or wrong type
    }


    //==============================================================================
    // Whenever the APVTS notifies the listener that a parameter has received a new
    // value, all this does is set the "parametersChanged" boolean to true. This is
    // thread-safe, since it's an atomic variable.
    std::atomic<bool> parametersChanged {false};
    juce::String lastChangedParamID {};
    void valueTreePropertyChanged(juce::ValueTree& tree, const juce::Identifier&) override {
        // DBG("Detected Parameter Change.");
        parametersChanged.store(true);
        lastChangedParamID = tree.getProperty("id").toString();
    }


    //==============================================================================
    // This "updateParameters()" method is where you'll do all the necessary
    // calculations using the new parameter values.
    void updateParameters();
    void refreshParameters();
    void setInitalConditions();


    //==============================================================================
    // "apvts.getRawParameterValue(...)->load()" performs a relatively inefficient
    // parameter lookup by comparing string, since parameter identifiers are
    // "juce::String" objects under the hood. You'll be doing this lookup hundreds
    // of times per second, since it happens at the start of "processBlock()".
    // It's more efficient use the AudioParameter object directly. The following
    // variables are used to catch the parameters coming from apvts. Each of them
    // are one by one connected to the parameters we defined at the beginning
    // namespace "MyParameterID".
    juce::AudioParameterFloat* inputGainParameter {};
    float inputGain { 0.0f };
    float inputAmp { 1.0f };

    juce::AudioParameterFloat* linParameter {};
    float alpha_lin { 0.5f };

    juce::AudioParameterFloat* volParameter {};
    float alpha_vol { 0.5f };


    //==============================================================================
    float Fs {};
    float T {};
    void set_fs(double sampleRate);


    //==============================================================================
    // Circuit Parameters
    const float C71  = 4.7e-6f;
    float C71_z {};

    const float Rin   = 22e3f;
    const float R71_z = Rin;

    const float R72   = 50e3f;
    const float R72_z = R72;

    const float C81  = 4.7e-6f;
    float C81_z {};

    const float R81   = 22e3f;
    const float R81_z = R81;

    const float Rout  = 100e3f;
    const float R91_z = Rout;

    const float C91  = 4.7e-6f;
    float C91_z {};

    const float Rvol = 50e3f;
    float R92_z {};
    float R93_z {};

    const float C94  = 4.7e-6f;
    float C94_z {};


    //==============================================================================
    // Circuit State Variables (Kirchhoff Domain)
    const float Vcc = { 9.0f };
    const float C71_ic = { -0.820414245128631591796875f };        // @ lin: 50%, vol: 50%
    const float C81_ic = { 0.17600531876087188720703125f };       // @ lin: 50%, vol: 50%
    const float C91_ic = { 0.0010983161628246307373046875f };     // @ lin: 50%, vol: 50%
    const float C94_ic = { 8.82324504852294921875f };             // @ lin: 50%, vol: 50%

    // Circuit State Variables (Wave Domain)
    std::array<float, 2> C71_a {};
    std::array<float, 2> C71_b {};
    std::array<float, 2> C81_a {};
    std::array<float, 2> C81_b {};
    std::array<float, 2> C91_a {};
    std::array<float, 2> C91_b {};
    std::array<float, 2> C94_a {};
    std::array<float, 2> C94_b {};


    //==============================================================================
    Eigen::Vector<float, 3> S71_a {};
    Eigen::Vector<float, 3> S71_b {};
    Eigen::Vector<float, 3> P72_a {};
    Eigen::Vector<float, 3> P72_b {};
    Eigen::Vector<float, 3> P81_a {};
    Eigen::Vector<float, 3> P81_b {};
    Eigen::Vector<float, 3> S91_a {};
    Eigen::Vector<float, 3> S91_b {};
    Eigen::Vector<float, 3> P92_a {};
    Eigen::Vector<float, 3> P92_b {};
    Eigen::Vector<float, 3> S93_a {};
    Eigen::Vector<float, 3> S93_b {};
    Eigen::Vector<float, 3> S94_a {};
    Eigen::Vector<float, 3> S94_b {};

    Eigen::Vector<float, 3> S71_Zv {};
    Eigen::Vector<float, 3> P72_Zv {};
    Eigen::Vector<float, 3> P81_Zv {};
    Eigen::Vector<float, 3> S91_Zv {};
    Eigen::Vector<float, 3> P92_Zv {};
    Eigen::Vector<float, 3> S93_Zv {};
    Eigen::Vector<float, 3> S94_Zv {};

    const Eigen::Vector<float, 3> S71_Bv {  1.0f, 1.0f, -1.0f };
    const Eigen::Vector<float, 3> P72_Qv {  1.0f, 1.0f,  1.0f };
    const Eigen::Vector<float, 3> P81_Qv {  1.0f, 1.0f,  1.0f };
    const Eigen::Vector<float, 3> S91_Bv {  1.0f, 1.0f,  1.0f };
    const Eigen::Vector<float, 3> P92_Qv {  1.0f, 1.0f,  1.0f };
    const Eigen::Vector<float, 3> S93_Bv {  1.0f, 1.0f,  1.0f };
    const Eigen::Vector<float, 3> S94_Bv { -1.0f, 1.0f,  1.0f };

    Eigen::Matrix<float, 3, 3> S71_S {};
    Eigen::Matrix<float, 3, 3> P72_S {};
    Eigen::Matrix<float, 3, 3> P81_S {};
    Eigen::Matrix<float, 3, 3> S91_S {};
    Eigen::Matrix<float, 3, 3> P92_S {};
    Eigen::Matrix<float, 3, 3> S93_S {};
    Eigen::Matrix<float, 3, 3> S94_S {};

    void update_port7_scattering_matrix(float fs);
    void update_port8_scattering_matrix(float fs);
    void update_port9_scattering_matrix(float fs, float alpha_vol_);

    template<size_t NumPortsTotal>
    Eigen::Matrix<float, NumPortsTotal, NumPortsTotal>
    get_scattering_matrix_series(
            const Eigen::Vector<float, NumPortsTotal>& Bv,
            const Eigen::Vector<float, NumPortsTotal>& Zv) const;

    template<size_t NumPortsTotal>
    Eigen::Matrix<float, NumPortsTotal, NumPortsTotal>
    get_scattering_matrix_parallel(
            const Eigen::Vector<float, NumPortsTotal>& Qv,
            const Eigen::Vector<float, NumPortsTotal>& Zv) const;


    //==============================================================================
    Eigen::Vector<float, 4> RJ_af {};
    Eigen::Vector<float, 2> RJ_bf {};
    Eigen::Vector<float, 6> RJ_ab {};
    Eigen::Vector<float, 3> RJ_bb {};
    size_t lin_index { 50 };
    size_t vol_index { 50 };

    HighFive::File SZ_h5File;
    std::array<std::array<Eigen::Matrix<float, 2, 4>, 101>, 101> RJ_Sf {};
    std::array<std::array<Eigen::Matrix<float, 3, 6>, 101>, 101> RJ_Sb {};
    std::array<std::array<Eigen::Vector<float, 3>, 101>, 101> RJ_Zr {};
    Eigen::Matrix<float, 2, 4> RJ_Sf_current {};
    Eigen::Matrix<float, 3, 6> RJ_Sb_current {};
    Eigen::Vector<float, 3> RJ_Zr_current {};
    bool load_RJ_scattering_matrix(float fs);

    enum FS {
        RATE_44100,   // 对应 44.1 kHz
        RATE_48000,   // 对应 48 kHz
        RATE_88200,   // 对应 88.2 kHz
        RATE_96000,   // 对应 96 kHz
        RATE_176400,  // 对应 176.4 kHz
        RATE_192000   // 对应 192 kHz
    };

    const std::unordered_map<float, size_t> fs_dict = {
            {44.1e3,  FS::RATE_44100},
            {48e3,    FS::RATE_48000},
            {88.2e3,  FS::RATE_88200},
            {96e3,    FS::RATE_96000},
            {176.4e3, FS::RATE_176400},
            {192e3,   FS::RATE_192000}
    };


    //==============================================================================
    // BJT model parameters (Wave Domain)
    Eigen::Matrix<float, 2, 2> BJT_Z {};
    Eigen::Vector<float, 5> BJT_in {};
    Eigen::Vector<float, 2> BJT_out {};


    //==============================================================================
    // Neural Network Model
    void loadNN();

    template <typename BatchNormType>
    void loadBatchNorm1d(const HighFive::File& h5File, const std::string& layerName, BatchNormType& batchNorm1d);

    template <typename DenseType>
    void loadDense(const HighFive::File& h5File, const std::string& layerName, DenseType& dense);

    RTNeural::BatchNorm1DT<float, 5> gn1 {};
    RTNeural::DenseT<float, 5, 16> linear1 {};
    RTNeural::ReLuActivationT<float, 16> relu1 {};
    RTNeural::DenseT<float, 16, 16> linear2 {};
    RTNeural::ReLuActivationT<float, 16> relu2 {};
    RTNeural::DenseT<float, 16, 2> linear3 {};


    //==============================================================================



    //==============================================================================
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (AudioPluginAudioProcessor)
};

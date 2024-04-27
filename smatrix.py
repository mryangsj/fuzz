import numpy as np

def series(
    Z: np.ndarray | list | tuple
) -> np.ndarray:
  """
  Calculate the scattering matrix of a series junction based on the given free-parameter impedance vector Z. This function uses the method of general wave domain (WD) Juction to derive the scattering matrix.

  Parameters:
  Z (np.ndarray, list, tuple): 
    A free-parameter impedance vector. Each element of the vector represents the free-parameter impedance of a port.

  Returns:
  np.ndarray: The Scattering Matrix of the series junction.
  """

  # STEP 1: Formate the Z vector
  try:
    Z = np.array(Z, dtype=float).reshape(1, -1) # Attempt to convert Z to a 1xn numpy array
  except ValueError:
    raise ValueError("The free-parameter impedance vector Z cannot be converted to a 1xn numpy array.")
  if Z.size == 0:  # Check if Z is empty
    raise ValueError("The free-parameter impedance vector Z is empty.")
  if np.all(Z == 0):  # Check if all elements in Z are 0
    raise ValueError("The free-parameter impedance vector Z cannot be all zeros.")
  
  # STEP 2: Formulate fundamental loop matrix
  num_ports = Z.shape[1]  # Get the number of ports
  B = np.ones((1, num_ports))

  # STEP 3: Calculate the scattering matrix using the general WD junction method
  S = generalB(Z, B)
  
  # STEP 4: Return the scattering matrix
  return S


def parallel(
    Z: np.ndarray | list | tuple
) -> np.ndarray:
  """
  Calculate the scattering matrix of a parallel junction based on the given free-parameter impedance vector Z. This function uses the method of general wave domain (WD) Jucntion to derive the scattering matrix.
  
  Parameters:
  Z (np.ndarray, list, tuple): A free-parameter impedance vector. Each element of the vector represents the free-parameter impedance of a port.

  Returns:
  np.ndarray: The Scattering Matrix of the parallel junction.
  """

  # STEP 1: Formate the Z vector
  try:
    Z = np.array(Z, dtype=float).reshape(1, -1) # Attempt to convert Z to a 1xn numpy array
  except ValueError:
    raise ValueError("The free-parameter impedance vector Z cannot be converted to a 1xn numpy array.")
  if Z.size == 0:  # Check if Z is empty
    raise ValueError("The free-parameter impedance vector Z is empty.")
  if np.all(Z == 0):  # Check if all elements in Z are 0
    raise ValueError("The free-parameter impedance vector Z cannot be all zeros.")
  
  # STEP 2: Formulate fundamental cut-set matrix
  num_ports = Z.shape[1]  # Get the number of ports
  Q = np.ones((1, num_ports))  # Build the fundamental loop matrix of parallel junction
  
  # STEP 3: Calculate the scattering matrix using the general WD junction method
  S = generalQ(Z, Q)
  
  # STEP 4: Return the scattering matrix
  return S

def generalB(
    Z: np.ndarray | list | tuple,
    B: np.ndarray | list | tuple
) -> np.ndarray:
  """
  Calculate the scattering matrix of a general junction based on the given free-parameter impedance vector Z and the fundamental loop matrix B.
  
  Parameters:
  Z (np.ndarray, list, tuple): A free-parameter impedance vector. Each element of the vector represents the free-parameter impedance of a port.
  B (np.ndarray, list, tuple): The fundamental loop matrix of the junction.

  Returns:
  np.ndarray: The Scattering Matrix of the general junction.
  """

  # STEP 1: Formate the Z vector to a diagonal matrix
  try:
    Z = np.array(Z, dtype=float).reshape(1, -1) # Attempt to convert Z to a 1xn numpy array
  except ValueError:
    raise ValueError("The free-parameter impedance vector Z cannot be converted to a 1xn numpy array.")
  if Z.size == 0:  # Check if Z is empty
    raise ValueError("The free-parameter impedance vector Z is empty.")
  if np.all(Z == 0):  # Check if all elements in Z are 0
    raise ValueError("The free-parameter impedance vector Z cannot be all zeros.")
  
  Z = np.diagflat(Z)  # Convert Z to a diagonal matrix
  
  # STEP 2: Formate the B matrix
  try:
    B = np.array(B, dtype=float) # Attempt to convert B to a numpy array
  except ValueError:
    raise ValueError("The fundamental loop matrix B cannot be converted to a numpy array.")
  if B.size == 0:  # Check if B is empty
    raise ValueError("The fundamental loop matrix B is empty.")
  if np.all(B == 0):  # Check if all elements in B are 0
    raise ValueError("The fundamental loop matrix B cannot be all zeros.")
  
  # STEP 3: Calculate the scattering matrix using the general WD junction method
  num_ports = Z.shape[0]  # Get the number of ports
  S = np.eye(num_ports) - 2 * Z @ B.T @ np.linalg.inv(B @ Z @ B.T) @ B
  
  # STEP 4: Return the scattering matrix
  return S


def generalQ(
    Z: np.ndarray | list | tuple,
    Q: np.ndarray | list | tuple
) -> np.ndarray:
  """
  Calculate the scattering matrix of a general junction based on the given free-parameter impedance vector Z and the fundamental cut-set matrix Q.
  
  Parameters:
  Z (np.ndarray, list, tuple): A free-parameter impedance vector. Each element of the vector represents the free-parameter impedance of a port.
  Q (np.ndarray, list, tuple): The fundamental cut-set matrix of the junction.

  Returns:
  np.ndarray: The Scattering Matrix of the general junction.
  """

  # STEP 1: Formate the Z vector to a diagonal matrix
  try:
    Z = np.array(Z, dtype=float).reshape(1, -1) # Attempt to convert Z to a 1xn numpy array
  except ValueError:
    raise ValueError("The free-parameter impedance vector Z cannot be converted to a 1xn numpy array.")
  if Z.size == 0:  # Check if Z is empty
    raise ValueError("The free-parameter impedance vector Z is empty.")
  if np.all(Z == 0):  # Check if all elements in Z are 0
    raise ValueError("The free-parameter impedance vector Z cannot be all zeros.")
  
  Z = np.diagflat(Z)  # Convert Z to a diagonal matrix
  
  # STEP 2: Formate the Q matrix
  try:
    Q = np.array(Q, dtype=float) # Attempt to convert Q to a numpy array
  except ValueError:
    raise ValueError("The fundamental cut-set matrix Q cannot be converted to a numpy array.")
  if Q.size == 0:  # Check if Q is empty
    raise ValueError("The fundamental cut-set matrix Q is empty.")
  if np.all(Q == 0):  # Check if all elements in Q are 0
    raise ValueError("The fundamental cut-set matrix Q cannot be all zeros.")
  
  # STEP 3: Calculate the scattering matrix using the general WD junction method
  num_ports = Z.shape[0]  # Get the number of ports
  Z_inv = np.linalg.inv(Z)
  S = 2 * Q.T @ np.linalg.inv(Q @ Z_inv @ Q.T) @ Q @ Z_inv - np.eye(num_ports)
  
  # STEP 4: Return the scattering matrix
  return S
import numpy as np

c = 343

def compute_reverb_weighting(room, α):
    """
    @brief  computes the reverberation-weighting used in the paper.
    @param  the room object
    @return the weighting
    """
    A = 0
    for wall in room.walls:
        A = A + room.wall_area(wall)
    V = room.get_volume()
    D = 2.5
    Q = 1
    T = 24*np.log(10)*V/c/α/A
    σ = 16*np.pi*D**2*(A*T-0.163*V)/0.163/Q/A/V
    γ = 1/(1 + σ)
    return γ

"""
The Analysis of the GCC-PHAT
"""

def compute_tdoa(x_0, x_1, f_s, γ):
    """
    @brief  Computes the time-difference of arrival (TDOA) from two frames from two
            microphones by GCC-PHAT.
    @param  x_0 (F, 1): the first channel
    @param  x_1 (F, 1): the second channel
    @param  f_s: the sampling frequency in Hz
    @return float: the TDOA  
    """

    # Work out the frame-length.
    F = x_0.shape[0]

    # Compute the FFT of both signals.
    X_0 = np.fft.fft(x_0)
    X_1 = np.fft.fft(x_1)

    # This fixes a bug where the last bin is zero.
    if X_0[4096] == 0j:
        X_0[4096] = X_0[4095]
    if X_1[4096] == 0j:
        X_1[4096] = X_1[4095]

    # Compute the mean across the whole spectrum and make a weighting based 
    # thereon.
    X_m = (np.absolute(X_0) + np.absolute(X_1))/2
    X_n = np.mean(X_m)
    w = (X_m >= X_n) + (X_m < X_n)*(X_m/X_n)**0.3

    # Compute the spectral weighting, i.e. PHAT.
    φ = 1 / (γ * np.absolute(X_0) * np.absolute(X_1))
    # φ = w**2 / (γ * np.absolute(X_0) * np.absolute(X_1))
    # φ = 1

    # Compute the GCC.
    χ = np.multiply(φ, np.multiply(X_0, np.conj(X_1)))
    R_01 = np.fft.fftshift(np.real(np.fft.ifft(χ)))

    # Compute the TDOA and the azimuth.
    τ = np.linspace(-(F+1)/f_s/2, (F-1)/f_s/2, F)
    Δt_01 = τ[np.argmax(R_01)]

    return Δt_01

def estimate_position_from_one(r_0, Δt, i, r_m):
    """
    @brief  Estimates the position of the source from one set of equations given 
            a microphone as a reference.
    @param  r_0 (3, 1): the initial position
    @param  Δt (M, M): the computed TDOAs for all pairs
    @param  i: the index of the microphone not taken
    @return ndarray (3, 1): the estimated position
    """

    K = 100

    # Calculate the distances to the source for each pair of microphones. 
    Δd_0 = c*np.delete(Δt[:,0], [0, i]).T
    
    def f(r):
        """
        @brief  Computes the localisation-model.
        @param  r (3, 1): the estimated position of the source
        @return ndarray (M-1, 1): the error, i.e. the difference
        """
        d = np.linalg.norm(r - r_m.T, axis = 1)
        return np.delete(d, [0, i]) - d[0] - Δd_0
    
    def J(r):
        """
        @brief  Computes the Jacobian of the localisation-model.
        @param  r (3, 1): the estimated position of the source
        @return ndarray (M-1, 3): the Jacobian for the given position of the source
        """
        d = np.linalg.norm(r - r_m.T, axis = 1)
        return (
            (r - np.delete(r_m.T, [0, i], axis = 0)) 
            / np.delete(d, [0, i]).reshape((3,1))
            - (r - r_m.T[0]) / d[0]
        )

    # Iterate through the localisation-model to estimate the position of the 
    # source.
    for k in range(K):
        r_0 = r_0 - np.linalg.solve(J(r_0), f(r_0))
        # If the position is far off, then the iterative algorithm has failed.
        if np.linalg.norm(r_0) >= 10**2:
            raise OverflowError("Failure to compute iteration")
    
    return r_0

def estimate_position_from_all(x, f_s, r_m):
    """
    @brief  Estimates the position of the source from many sets of equations.
    @param  x (F, M): the frame with all channels
    @param  f_s: the sampling frequency in Hz
    @return ndarray (3, 1): the estimated position
    """

    M = x.shape[1]

    # Compute the TDOAs for all pairs.
    Δt = np.zeros((M, M))
    for i in range(M):
        for j in range(M):
            if j == i:
                continue
            Δt[i, j] = compute_tdoa(x[:,i], x[:,j], f_s)

    # Average the diagonal pairs, e.g. x[0,1] and x[1,0].
    Δt = (Δt - Δt.T)/2
    
    # Compute the initial position from the partitions.
    # Chen & Xu do not give much information about how they computed this. From 
    # their brief example, I have assumed that it checking whether the TDOAs on 
    # both sides of a chosen microphone on the pyramid's base are positive, 
    # i.e. Δt_10 > 0, Δt_30 > 0. If they are, then see which adjacent 
    # microphone has the shortest TDOA. The partition lies in the sector 
    # between the chosen microphone and the bisecting line between it and the 
    # adjacent microphone.
    r_0 = np.zeros((3,))
    for m in range(1, M):
        # If the TDOAs on both sides of the m-th microphone are non-zero, then
        # check the other neighbouring two.
        if Δt[(m - 0) % (M - 1) + 1, m] >= 0 \
        and Δt[(m - 2) % (M - 1) + 1, m] >= 0:
            # Choose one of the neighbouring microphones
            if Δt[(m - 0) % (M - 1) + 1, m] < Δt[(m - 2) % (M - 1) + 1, m]:
                l = (m - 0) % (M - 1) + 1
            else:
                l = (m - 2) % (M - 1) + 1
            # Find the bisecting line between the two microphones.
            r_mid = (r_m[:,m] + r_m[:,l]) / 2
            # Span outwards to find a good enough initial position between the 
            # two lines.
            r_0 = 2.0*(r_m[:,m] + r_mid)

    # Estimate the source's position as the mean of five estimates each with 
    # reference to each microphone.
    r = np.zeros((3, M))
    scaling = 1.0
    for m in range(1, M):
        # r[:,m] = estimate_position_from_one(r_0, Δt, m)
        try:
            r[:,m] = estimate_position_from_one(scaling*r_0, Δt, m, r_m)
        except OverflowError as error:
            if scaling == 3.0:
                raise OverflowError("Failure to compute iteration")
            scaling = scaling + 0.5
            m = m - 1
            continue
        scaling = 1.0
    r_F = np.sum(r, axis = 1) / (M-1)

    return r_F.reshape((3,1))

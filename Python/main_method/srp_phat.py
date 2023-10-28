import numpy as np
import matplotlib.pyplot as plt

import time

class srp_phat:
    def __init__(self, r_m, f_s, F, c = 343):
        """
        @brief  estimates the location by SRP-PHAT.
        @param  r_m (3, M): the coordinates of the microphones with respect to 
                the array,
        @param  f_s: the sampling-frequency in Hz,
        """
        self.r_m = r_m
        self.M = r_m.shape[1]
        self.f_s = f_s
        self.F = F
        self.c = c
        self.Y = np.zeros(F)
        self.Y_N = np.zeros(F)
        self.n_frames = 0
    
    def make_circular_grid(self):
        G = 144
        θ = np.deg2rad(2.5)
        v = np.zeros((3,G))
        v[:,0] = np.array([1,0,0])
        R = np.array([
            [np.cos(θ), -np.sin(θ), 0],
            [np.sin(θ), np.cos(θ), 0],
            [0, 0, 0]
        ])
        for g in range(1,G):
            v[:, g] = np.matmul(R, v[:,g-1])

        # fig = plt.figure()
        # ax = fig.add_subplot(projection='3d')
        # ax.scatter(v[0,:], v[1,:], v[2,:])
        # plt.show()

        return v

    def make_spherical_grid(self, N_ev):
        """
        @brief  Makes the spherical grid as described in the papers by the French 
                Canadians  namely by tesselating an icosahedron.
        @param  N_ev: the number of evolutions,
        @return ndarray (3,?): an array of the points,
        """

        # Make an icosahedron with which to being.
        ψ_r = (1 + np.sqrt(5))/2
        v = np.array([
            [0,     1,      ψ_r ],
            [1,     ψ_r,    0   ],
            [ψ_r,   0,      1   ],
            [0,     -1,     ψ_r ],
            [-1,    ψ_r,    0   ],
            [ψ_r,   0,      -1  ],
            [0,     1,      -ψ_r],
            [1,     -ψ_r,   0   ],
            [-ψ_r,  0,      1   ],
            [0,     -1,     -ψ_r],
            [-1,    -ψ_r,   0   ],
            [-ψ_r,  0,      -1  ]
        ]).T

        # Normalise the points.
        v = v / np.linalg.norm(v, axis = 0).reshape((1,12))

        # Calculate the shortest distance between two points.
        d_span = np.linalg.norm(v[:,0] - v[:,3])

        # Tesselate the grid for however many evolutions.
        for k in range(N_ev):
            # For each point and for each of its nearest neighbours and make 
            # another point between them.
            n_v = v.shape[1]
            for i in range(n_v):
                for j in range(i+1,n_v):
                    # See if the two points are nearest neighbours by seeing 
                    # whether they are withing the shortest distance.
                    d = np.linalg.norm(v[:,i] - v[:,j])
                    if d >= 1.5*d_span/(2**k):
                        continue
                    # Make the new point and add it to the array.
                    u = (v[:,i] + v[:,j]).reshape((3,1))
                    u = u / np.linalg.norm(u)
                    v = np.hstack((v,u))

        # print("{} points".format(v.shape[1]))

        # fig = plt.figure()
        # ax = fig.add_subplot(projection='3d')
        # ax.scatter(v[0,:], v[1,:], v[2,:])
        # plt.show()

        return v

    def compute_tdoa(self, u):
        """
        @brief  Computes the TDOA for every pair given the direction.
        @param  u (3,1): the direction,
        @return ndarray (M,M): the TDOAs in number of samples,
        """
        τ = np.zeros((self.M,self.M))

        # For each pair of microphones, compute the TDOA given by the formula in 
        # Valin et al. 2007.
        for i in range(self.M):
            for j in range(self.M):
                if i == j:
                    continue
                τ[i,j] = np.dot(self.r_m[:,j]-self.r_m[:,i], u) \
                    * self.f_s / self.c
        return τ

    def compute_tdoa_grid(self, grid):
        """
        @brief  Computes the TDOA for every pair and for every direction in the 
                grid.
        @param  grid (3,G): the sperical grid,
        @return ndarray (G,M,M): the TDOAs in number of samples,
        """

        G = grid.shape[1]
        τ = np.zeros((G,self.M,self.M))

        # For every direction in the grid, compute the TDOAs.
        for g in range(G):
            τ[g] = self.compute_tdoa(grid[:,g])

        return τ

    def compute_gcc_phat(self, X_0, X_1):
        """
        @brief  Computes the GCC-PHAT from two frames from two microphones.
        @param  X_0 (F, 1): the first channel
        @param  X_1 (F, 1): the second channel
        @return ndarray (F, 1): the GCC-PHAT
        """

        # Work out the frame-length.
        F = X_0.shape[0]

        # Make a weighting based the spectral noise.
        # w = (self.Y <= self.Y_N) + (self.Y > self.Y_N)*(self.Y/self.Y_N)**0.99
        w = 1

        # Compute the spectral weighting, i.e. PHAT.
        # ψ = 1 / (np.absolute(X_0) * np.absolute(X_1))
        ψ = w**2 / (1 * np.absolute(X_0) * np.absolute(X_1))
        # ψ = 1 / (np.absolute(X_0) * np.absolute(X_1))
        # ψ = ψ * (np.linspace(0, self.f_s/2, F) < 10000)
        # ψ = 1

        # Compute the GCC.
        χ = np.multiply(ψ, np.multiply(X_0, np.conj(X_1)))
        R_01 = np.real(np.fft.ifft(χ))

        return R_01

    def compute_beam_energy(self, x, g, R, τ):
        """
        @brief  Computes the beamformer's energy along the g-th direction in the 
                grid.
        @param  x (F, M): the frame with all channels,
        @param  g: the index of the direction along which to compute,
        @param  R (M, M, F): the array of all GCC-PHATs,
        @param  τ (G, M, M): the array of all TDOAs for all directions,
        @return float: the energy,
        """
        
        # E = 0
        # for i in range(self.M):
        #     for j in range(self.M):
        #         if i == j:
        #             continue
        #         E = E + R[i,j,τ[g,i,j]]
        τ_space = τ[g,:,:].reshape((1,self.M**2))
        R_space = R.reshape((self.M**2,R.shape[2]))
        E = np.sum(R_space[np.arange(self.M**2),τ_space])
        return E
    
    def calculate_noise(self, x):
        """
        @brief  Calculates the spectral noise as an average in time of the mean 
                spectral density across all microphones.
        @param  x (F, M): the frame with all channels,
        """
        # Compute the FFT of all signals.
        X = np.fft.fft(x, axis=0)

        # Compute the mean spectral density across all microphones.
        self.Y = np.mean(np.abs(X), axis=1)

        # Calculate the time-average with a falling horizon to the beginning of 
        # time.
         
        # The French Canadians did not go into depth about what they 
        # exactly meant by the time-average. Here, it is taken as the total 
        # average rather than a moving average.
        self.Y_N = (self.Y_N*self.n_frames + self.Y)/(self.n_frames+1)
        self.n_frames = self.n_frames + 1

    def estimate_direction(self, x, τ, grid, map = False, ax = None):
        """
        @brief  Estimates the position of the source from one set of equations given 
                a microphone as a reference.
        @param  x (F, M): the frame with all channels,
        @param  τ (G, M, M): the array of all TDOAs for all directions,
        @param  grid (3, G): the array of directions on the grid,
        @param  map: whether a map is wanted,
        @param  ax: the axis for a map,
        @return tuple: the estimated azimuth and elevation in radians,
        """

        # self.calculate_noise(x)

        # Compute the FFT of all signals.
        # chrono = time.time()
        X = np.fft.fft(x, axis=0)
        # print("FFT done in ", time.time() - chrono, " s")

        # Compute the GCC-PHATs for all pairs.
        F = x.shape[0]
        R = np.zeros((self.M,self.M,F))
        # chrono = time.time()
        for i in range(self.M):
            for j in range(self.M):
                if i == j:
                    continue
                R[i,j,:] = self.compute_gcc_phat(X[:,i], X[:,j])
        # print("GCC done in ", time.time() - chrono, " s")
        τ_int = (np.rint(τ)).astype(int)

        # For each direction in the grid, compute the energy and see whether it is 
        # the largest energy.
        G = grid.shape[1]
        g_max = 0
        E_max = 0
        # chrono = time.time()
        for g in range(G):
            E = self.compute_beam_energy(x, g, R, τ_int)
            if E >= E_max:
                g_max = g
                E_max = E
        # print("Beam done in ", time.time() - chrono, " s")

        u = grid[:, g_max]
        θ = np.arctan2(u[1], u[0])
        ψ = np.arctan2(np.linalg.norm(u[0:2]), u[2])

        # If a map is wanted, then compute the energy-map again (this is to 
        # keep this code away from the rest which needs to be quick), build a 
        # polar coordinate-space, and for each coordinate, find the nearest 
        # direction in the grid and sample its energy. Then, plot the heat-map.
        if map == True:
            # Compute the energy-map again.
            E_grid = np.zeros(G)
            for g in range(G):
                E_grid[g] = self.compute_beam_energy(x, g, R, τ_int)
            # Make a linear space of polar coordinates.
            θ_length = 200*2
            ψ_length = 100*2
            E_space =np.zeros((θ_length,ψ_length))
            θ_lin = np.linspace(0, np.pi, θ_length)
            ψ_lin = np.linspace(-np.pi, np.pi, ψ_length)
            θ_space, ψ_space = np.meshgrid(θ_lin, ψ_lin)
            # For each coordinate, find the nearest direction in the grid and 
            # takes its corresponding energy to plot.
            for i_θ in range(θ_length):
                for i_ψ in range(ψ_length):
                    # Convert the polar coordinates to cartesian.
                    u_pixel = np.array([
                        [np.sin(θ_lin[i_θ])*np.cos(ψ_lin[i_ψ])],
                        [np.sin(θ_lin[i_θ])*np.sin(ψ_lin[i_ψ])],
                        [np.cos(θ_lin[i_θ])]
                    ])
                    # Find the difference.
                    v_d = np.linalg.norm(u_pixel-grid, axis=0)
                    # Find the index of the nearest direction in the grid.
                    g_pixel = np.argmin(v_d)
                    E_space[i_θ,i_ψ] = E_grid[g_pixel]
            # Plot the grid of energies.
            if ax == None:
                plt.pcolormesh(
                    np.rad2deg(ψ_space), 
                    np.rad2deg(θ_space), 
                    E_space.T
                )
                plt.colorbar()
            else:
                ax.pcolormesh(
                    np.rad2deg(ψ_space), 
                    np.rad2deg(θ_space), 
                    E_space.T
                )

        return (θ, ψ)
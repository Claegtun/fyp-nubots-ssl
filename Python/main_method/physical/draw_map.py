import matplotlib.pyplot as plt
import numpy as np

import sys
sys.path.append("./main_method")

import srp_phat

# Set the font.
plt.rcParams.update({'font.sans-serif':'FreeSerif'})

# Build the cube of microphones 15.5 cm wide.
centre = np.array([[5, 5, 1]]).T
r_m = np.array([
    [1,1,1],    # NE top
    [1,1,-1],   # NE bottom
    [-1,1,1],   # NW top
    [-1,1,-1],  # NW bottom
    [1,-1,1],   # SE top
    [1,-1,-1],  # SE bottom
    [-1,-1,1],  # SW top
    [-1,-1,-1]  # SW bottom
]).T * 91.2*10**(-3) / 2

f_s = 24*10**3

"""
The Beamforming Method
"""

# Get the number of samples in the frame.
F = 1024

# Make the object for the SRP-PHAT.
estimator = srp_phat.srp_phat(r_m, f_s, F)
# Compute the grid and all TDOAs.
grid = estimator.make_spherical_grid(4)

fig, ax = plt.subplots()

E_grid = np.array([
    0,1.271779,1.504823,2.088671,2.252486,1.046939,2.172803,2.039944,0.818169,1.859768,2.093136,1.162512,0.784933,0.638856,1.073744,1.691156,0.093300,0.414775,3.228866,1.249099,0.991409,0.779131,2.021278,1.246109,0.802521,1.034414,0.800780,0.869583,0.890559,1.018566,1.476518,0.706822,0.560642,0.817024,1.689991,1.101105,1.801025,0.893215,1.875962,1.145646,0.936479,0.781631,0.931349,0.779438,0.780702,0.245399,0.520791,0.406362,1.299293,2.361329,1.682283,0.739401,-0.054336,0.259139,0.443026,0.516838,0.648102,0.462568,0.481672,0.452482,0.536903,0.206787,-0.127786,1.701014,0.928205,1.466473,1.627025,-0.136363,0.702133,1.022282,0.190086,-0.439457,0.190793,1.040174,1.362030,0.016846,0.307008,2.064815,-0.222407,1.092811,0.245238,1.813715,0.064285,1.079232,0.279954,0.657112,0.373407,0.260011,0.581716,-0.233166,0.470653,0.859330,1.034528,1.769177,1.888310,0.235016,1.234069,1.738826,0.379982,0.808359,-0.066412,-0.091558,0.148302,0.386156,0.669632,0.023823,0.359283,0.206073,1.391566,1.338933,0.580810,0.379616,-0.097751,0.326361,0.860080,1.097273,-0.023131,3.659287,1.117122,0.001651,0.180720,0.593169,1.166528,-0.455312,0.362952,0.656537,0.563565,-0.198587,0.869600,0.734662,0.427149,3.471105,0.374148,-0.004452,0.086744,0.852327,0.113190,0.563032,0.768369,-0.260068,3.643737,0.709802,0.783506,1.146911,1.041725,1.732064,-0.197063,0.740408,1.189968,-0.053291,0.916973,0.712170,0.867574,0.179641,0.294148,0.247053,1.305561,3.139516,0.851243,0.083703,0.714250,0.759942,1.546525,0.894539,0.719471,0.656094,1.545421,1.222493,0.793304,1.820446,1.876245,0.983827,0.707707,0.187244,0.713162,-0.206697,0.633943,1.032352,1.965883,0.995580,1.064984,0.446147,1.394042,1.398611,1.181223,1.309515,1.266475,0.933476,-0.018282,1.203650,0.839448,0.488471,1.590857,1.241858,1.120550,0.974220,1.294208,0.837270,1.869173,0.866329,1.167909,1.587428,1.110286,0.689215,0.606434,0.636907,0.018075,1.144213,1.146972,0.524571,0.811028,1.659043,1.276010,1.710575,1.117509,1.009328,1.024270,0.361757,0.735620,0.078095,0.409389,1.025571,0.388541,0.510110,0.163291,1.058409,0.233821,-0.042978,0.994335,0.614438,0.781110,0.129351,0.968481,0.368366,1.463040,0.937028,0.542017,1.134382,0.913508,0.265045,0.806356,0.863009,0.507333,-0.048929,-0.105615,1.099573,0.453690,0.570492,0.047229,0.715903,0.986908,0.095493,1.445514,1.111875,0.976679,1.015393,2.162313,1.040820,2.458763,1.682490,1.531446,1.038026,0.400321,0.806164,1.554327,1.524338,1.281252,1.692069,0.472193,1.047073,1.764883,1.617067,0.852249,0.751924,0.362489,0.275930,1.231876,1.395855,1.008097,0.441760,0.342377,0.886253,0.549013,0.437945,0.401840,0.970122,1.311433,1.013000,1.785375,2.276623,1.513434,0.064244,1.898617,2.593051,1.004083,0.310062,0.321003,1.176301,0.698724,0.677271,0.512538,0.327299,0.255351,0.605486,0.662059,0.902699,0.716681,1.442858,0.346498,0.014968,0.405899,0.469431,0.120450,0.715744,0.551871,0.993996,-0.074415,0.690488,1.315107,0.823118,1.160313,1.700097,1.954820,0.020964,1.992989,1.702984,1.437469,3.018363,1.376007,-0.202293,0.600867,0.494776,0.841699,1.367741,0.792149,0.694468,0.631384,0.870337,-0.257563,2.329161,1.111053,0.524591,0.560442,0.197505,0.489697,0.622257,0.195159,0.715264,0.867353,0.965890,0.939005,1.023901,0.691100,0.556928,0.601090,1.189500,0.519974,0.663972,0.831409,1.222123,0.535183,0.368997,0.735839,0.115040,0.263505,0.147506,1.627152,2.112399,0.876518,0.570633,1.713030,0.658042,1.367480,1.496195,1.184569,0.909377,0.535952,1.707975,0.955866,0.896385,1.069635,1.165950,1.171402,0.992198,0.958666,0.045897,0.412194,0.119256,0.345274,0.491512,0.588138,0.331651,0.307316,0.428867,1.526120,0.355324,1.334142,1.122106,2.249674,0.155117,1.039367,1.479791,0.801813,1.749515,0.114340,0.096622,0.774242,0.905258,0.411404,0.715343,1.607378,1.013500,0.850640,0.335089,0.249223,0.323924,0.983228,1.240397,0.341908,0.446844,0.452390,2.460332,0.436306,2.005832,1.239444,1.102543,1.323312,2.460502,2.269324,0.932632,0.512067,0.189472,0.471638,1.151985,0.361785,0.047988,-0.220563,1.044049,0.846270,1.554509,1.146264,0.092963,0.436621,1.423027,1.946790,1.259476,0.298903,0.767695,0.153651,0.913363,1.147085,1.368842,0.802912,0.570834,1.791888,0.719852,-0.051976,0.772918,0.486330,0.927144,0.191482,0.133033,0.684194,2.183069,0.558600,-0.215183,0.417129,1.015039,0.597688,1.797820,0.038505,0.253264,2.472173,0.127409,2.468415,0.925220,2.119707,0.786844,1.672749,0.061230,1.927303,0.420131,1.527118,1.111890,1.172510,0.538226,0.946775,2.209507,-0.308732,1.237053,0.638444,0.575172,1.568355,1.282447,0.058276,0.047700,1.913270,0.651972,0.556244,0.760593,0.469466,0.666073,1.044134,0.645297,0.792235,0.634897,0.822386,1.776521,2.774968,1.114678,1.815347,0.845433,0.745421,0.054598,1.238989,2.901645,0.367018,0.880569,-0.216549,0.621263,0.183022,0.749066,1.853428,1.004523,0.901816,0.417188,1.556041,1.260304,1.239942,0.570822,0.893801,-0.012819,0.475113,1.343432,0.100108,1.398765,2.049648,1.814726,0.818014,0.185448,0.790641,1.142009,0.743594,1.225706,1.740903,0.122210,-0.075238,0.669927,0.316487,0.064405,0.756330,0.584478,0.517704,2.268583,0.680007,0.466983,0.897549,0.910313,1.126552,3.022306,0.677425,2.142368,1.120741,0.015833,-0.260108,1.949012,0.294372,1.614011,0.990165,1.566826,0.355664,0.803970,-0.124260,1.035806,0.555025,1.612696,1.206492,0.796764,0.825350,0.220738,0.689819,1.674943,-0.278475,0.612228,0.591598,0.141529,0.718381,0.435240,0.399843,0.362768,0.951069,0.534688,0.671108,0.423844,0.628493,0.059122,0.320761,2.274973,1.396967,2.425609,-0.055534,0.266861,0.348727,0.763992,1.286851,0.177219,0.373972,0.239428,-0.377406,2.009762,1.843817,1.140796,0.440111,-0.219129,0.077010,0.643182,0.741155,0.539074,-0.130789,-0.202194,0.459511,1.107561,1.448215,1.471767,0.634080,0.494325,0.600637,0.320549,0.081631,0.320865,0.871233,1.038108,0.345486,0.174115,0.104145,2.110412,1.306933,2.327093,1.834655,-0.604082,0.145010,0.880807,2.373762,1.947866,0.630217,1.649738,1.043159,1.976737,2.182848,1.814955,2.127226,1.892085,0.630027,0.734889,1.533465,0.051903,1.498199,0.951809,1.729427,2.064767,1.042903,1.875465,1.810725,2.076151,2.209500,2.509913,2.516665,1.416448,0.971867,1.469187,1.111410,1.363241,1.641380,2.174134,1.166606,0.878824,2.365513,2.936220,2.234019,2.698014,2.067848,1.702235,1.107223,1.259823,0.779533,1.238774,-0.094050,1.121663,0.515684,1.432166,1.999106,1.561956,1.846621,1.295875,1.729061,2.356130,2.266180,1.029149,1.112211,0.145506,0.815214,0.541181,0.557959,2.051657,0.846718,0.797388,0.829527,0.341462,1.024948,0.380984,0.565177,0.043020,0.607041,0.098831,0.761488,1.350214,0.929800,0.477093,-0.194797,0.416332,1.057134,1.280683,1.608009,0.360063,1.232868,0.363741,0.188912,0.777352,0.509293,0.237263,0.650676,0.072411,1.338750,0.450580,1.141055,0.850645,1.908912,0.430510,2.096172,1.823744,1.658038,2.252999,1.797212,1.514001,1.724174,1.552340,2.071496,1.289010,2.350721,1.540860,0.898371,1.235357,0.925041,0.547308,1.707856,0.672881,1.619237,1.167414,0.842685,0.952462,0.749701,0.591584,1.348971,2.262219,1.253681,1.455641,0.876890,2.104514,1.418731,1.333007,0.703291,1.117506,0.961489,2.413282,1.791394,1.095767,0.392198,0.572991,0.946955,0.740843,0.717714,0.609588,0.723289,0.679704,0.730079,0.975096,1.172272,1.244242,0.627838,0.599409,1.398096,0.620249,0.773944,1.003055,1.412385,1.024902,0.420919,0.220745,1.285725,0.750329,0.768942,0.428718,2.221998,0.511737,1.650007,1.476148,0.750147,2.250277,1.058503,1.557711,0.723808,0.773782,1.069442,0.756813,0.740039,0.572579,0.957450,0.831277,0.962717,1.739752,1.100628,0.926000,1.069110,0.298231,0.153953,0.157937,0.318463,0.046492,0.270179,0.961225,0.620423,0.678906,0.585368,0.561495,0.691782,0.744458,1.392211,0.557157,0.881168,1.087839,0.943558,0.536429,1.048092,0.114913,1.326692,1.169096,1.862575,2.026813,1.506816,2.017430,2.394991,1.717176,1.012882,1.174627,1.260511,1.293707,1.582564,1.891215,1.926310,1.776128,1.590103,1.425310,1.562272,1.380512,0.615562,0.571711,1.049327,-0.274618,0.645463,0.758591,0.471961,0.026653,0.444933,0.052614,0.559885,0.457779,1.314940,0.634648,1.662515,1.214773,1.040597,0.706994,0.867383,1.041154,1.317116,0.755903,1.274347,0.755880,0.247273,0.783852,0.685869,0.308697,1.196013,0.179498,0.977149,1.421540,0.679324,1.136140,0.679214,1.457601,1.224024,0.319444,1.847166,0.848818,0.881801,0.212007,0.988800,0.488355,0.761338,0.973345,0.049706,0.779032,1.702931,0.823346,0.368471,1.273492,0.103142,0.062996,1.048448,1.709356,1.296861,1.637776,1.281647,2.840723,2.466176,2.059534,2.121746,2.777427,2.409205,2.064559,1.709223,1.950773,2.186990,1.646263,2.711576,1.112958,1.430837,0.990326,1.226712,1.166757,1.422380,0.455915,0.520758,0.544549,0.444128,0.566399,0.364356,0.397078,0.173822,1.422771,0.960970,-0.078195,1.518788,-0.000817,1.022260,1.152649,1.652713,1.557731,1.184622,0.233654,0.354106,0.986790,0.396430,0.440211,0.101166,0.986928,0.543227,1.228451,1.621901,1.077368,1.762095,0.246876,0.414654,0.401433,0.337891,0.399918,0.526083,0.383639,0.544969,0.937359,0.389030,0.964511,1.146439,0.793575,2.348804,1.231389,1.555662,1.440257,0.792063,0.821243,0.194564,0.298603,0.540751,0.567430,0.217646,0.991757,0.555878,0.897293,1.504667,0.506189,0.594988,0.323028,0.300052,0.139319,0.280255,0.661434,0.477140,0.284682,1.764398,1.102364,1.549524,1.709552,1.680652,1.914532,1.724864,0.207490,1.368848,0.861302,0.410830,-0.008417,1.544570,1.528102,1.422973,1.516096,0.628643,3.217375,1.837112,2.782962,1.454208,1.905218,1.656358,2.400810,-0.370791,1.560666,0.096848,1.704840,0.236902,1.642391,0.211745,1.498012,1.367179,0.449787,0.748813,1.490566,1.401276,1.264287,1.421909,1.950019,0.881162,0.940988,0.484389,0.817813,0.269620,1.488680,0.020636,1.009142,0.551413,0.403400,1.405985,0.867899,0.902250,0.508960,0.837023,1.341714,1.277130,0.720197,0.758762,0.411787,1.585710,0.215299,1.126421,1.421554,1.143878,0.731770,-0.015038,0.223385,0.670041,0.341779,0.961085,0.454253,0.284608,0.768058,0.503791,0.296526,0.283034,-0.212215,1.047988,0.729101,0.936271,0.604614,0.708956,0.532323,1.554368,2.820547,1.923347,1.953175,1.567388,2.888414,0.975739,0.044280,0.944395,0.982852,0.591335,0.473933,1.708060,1.216202,1.536614,1.337897,2.554381,0.566789,0.288989,0.242819,0.855604,0.739997,0.509870,-0.200524,1.752674,1.612647,1.594350,2.269214,1.609930,2.013906,0.594706,0.902458,0.879303,0.716257,1.303636,0.196357,1.700228,1.344061,1.475557,1.530223,1.303216,0.914115,0.354443,0.744159,2.214487,0.550509,1.108176,0.362465,0.384285,1.069670,1.301654,0.081383,1.651751,-0.098763,0.329640,1.437367,0.863740,0.432413,0.741465,0.728653,0.466080,0.780242,0.817310,0.772715,0.455849,0.615360,0.639961,1.243495,0.369690,0.808882,0.699057,0.726893,0.667409,0.329621,0.630358,1.345930,-0.136048,0.711762,1.291022,0.462119,1.234350,1.332773,0.456559,0.495868,0.590881,0.317667,0.653403,0.676040,0.202775,0.975725,1.014065,1.044567,1.013192,0.582972,1.103592,0.938576,2.150192,2.168998,2.103653,2.528659,1.943359,2.140584,1.538948,3.257611,1.665075,1.568631,1.531261,3.234716,1.234984,0.265874,1.055965,0.235594,0.695524,0.310717,1.017442,1.931790,1.635162,1.705805,1.993092,1.356581,0.487603,1.793949,1.711159,0.624480,1.687092,0.542373,1.445720,1.054001,2.234208,1.137135,0.498787,1.105219,1.408232,1.046063,0.681823,0.445793,0.845483,0.149085,0.869995,0.649011,0.609481,0.379035,0.799645,-0.081563,0.477913,1.306069,0.052882,1.071332,0.363643,1.666779,0.484069,0.893482,-0.031348,1.087454,0.033910,-0.166185,1.183593,1.232933,0.871632,0.463545,0.484080,-0.214170,0.237098,0.398413,0.433135,1.384163,0.102364,0.207628,0.356521,1.072855,1.024955,1.112550,-0.093705,0.078439,0.708320,1.028692,0.322863,1.122525,1.062549,0.798492,0.387422,1.130834,-0.011429,0.194992,0.855737,0.796207,0.745425,0.861290,1.509926,1.373777,0.790849,0.761805,1.102673,0.496622,1.033705,0.840475,0.830834,1.113186,0.650367,0.798741,1.044754,0.809824,0.694405,0.752833,0.451589,0.706851,1.082251,0.475017,0.423417,0.677937,0.052966,0.543669,0.273262,0.253824,0.139458,-0.261699,-0.090675,1.619002,0.620609,1.209860,0.372693,0.455347,0.405610,0.506895,-0.028403,0.314458,0.393446,0.573939,1.188423,0.998817,1.009897,1.465280,0.482655,0.711077,0.115997,1.259765,0.110392,0.697269,0.114939,0.467748,2.493120,1.673035,3.395746,3.635842,2.407009,2.044982,1.700752,1.534618,0.738050,1.156737,1.330828,0.962020,1.181856,0.067133,1.237007,0.121047,0.078423,0.086665,0.680175,0.202115,0.851184,0.185109,0.229972,0.212320,0.516111,1.203195,0.211290,0.807068,0.404527,1.141065,0.911834,2.461408,0.056129,1.399594,0.752848,1.417938,1.119160,0.425284,0.935614,0.391371,0.670093,0.642004,0.538454,0.510651,0.141208,-0.008558,1.016931,-0.007080,0.025372,0.688022,0.438178,0.958511,0.422240,0.110532,0.529059,1.544964,0.360350,0.721151,0.309977,0.428075,0.460692,-0.041761,0.120405,0.450209,0.138363,0.168496,1.226896,0.727566,2.256206,0.587109,1.224124,0.943079,1.123141,0.625059,0.994382,0.524464,0.771359,0.804335,0.564995,0.747044,1.652517,0.928371,0.283215,0.429110,2.760847,2.235325,4.025347,3.709260,2.585173,2.217486,0.262867,0.480659,1.361037,0.353705,0.774129,-0.222369,0.095374,0.509109,0.184462,0.978179,-0.056713,-0.195546,0.977819,0.250112,1.099293,-0.083065,0.420177,0.211073,0.631101,1.053479,0.172829,0.627694,0.424272,0.074707,0.102924,1.032112,0.384595,0.924790,0.181287,0.745691,0.770855,0.488360,0.582995,0.786866,1.440073,0.366455,0.412308,0.324022,0.751497,0.874943,0.861847,0.154077,0.152236,0.217027,0.532034,0.366133,0.043912,-0.349190,1.695647,3.120936,3.836556,3.645286,1.923505,2.694291,0.850154,1.341259,1.062583,0.882221,0.724821,0.964160,2.222744,0.660552,1.108941,0.681692,0.524582,0.691547,0.945433,0.940576,0.773397,1.256160,0.766355,0.253150,0.953688,1.162055,1.478734,1.036298,0.644401,0.808755,0.733375,1.071085,0.624367,0.365681,1.092020,0.516565,1.328770,0.426689,1.181390,0.861986,0.235274,0.260534,0.277786,0.523013,0.230625,-0.092058,0.207949,0.644613,1.225725,1.262908,0.533580,1.158902,1.761551,0.509150,-0.039216,0.038604,0.569779,-0.070237,0.083462,0.159132,0.635095,0.924259,1.039025,1.476052,0.520720,0.647497,1.206805,0.889758,1.451398,0.634721,0.660054,0.929764,0.717582,0.203458,0.980020,1.317875,0.402210,0.887519,0.001955,1.839228,0.437545,1.299726,0.197538,0.601354,0.066482,0.448618,0.587862,0.169910,0.074405,0.124493,1.631012,-0.076799,1.340440,-0.626206,0.313591,-0.322686,0.745034,1.473059,0.475557,1.267187,0.833468,1.509432,2.054923,2.501970,3.407486,3.185565,3.238737,2.764328,1.310322,1.389812,1.462206,0.635768,0.636803,1.063406,0.075366,0.738233,-0.199882,1.128649,-0.114481,0.120690,1.028291,0.525770,0.750643,0.497518,-0.094635,0.219869,0.461581,1.243464,0.187306,0.666245,0.354063,0.618528,1.379256,1.484100,1.293544,1.837918,1.319445,1.264747,0.713327,1.135793,1.511060,0.590169,1.425532,1.509616,0.929493,0.262029,0.915081,1.925757,2.075409,0.571034,1.426704,2.491507,1.568469,2.503215,1.940692,2.178590,2.048136,1.988575,2.583706,1.424222,1.333726,0.720884,0.569896,1.261544,0.490369,1.271035,0.443787,1.002642,0.568931,1.230748,0.895143,1.277279,-0.186499,0.157735,-0.087878,1.573131,0.568225,1.502031,0.641046,1.736390,0.266227,2.629581,1.216962,1.185697,1.639184,1.137791,2.096803,0.906908,0.206083,0.875528,1.343066,0.838380,2.239430,1.380610,0.744282,0.935688,2.766304,1.723740,2.014257,1.107924,1.044835,0.946345,3.076030,1.866107,2.300756,1.729277,2.873229,-0.229351,1.250451,-0.126107,1.210289,0.496949,-0.172316,0.575544,0.983153,1.557927,1.639231,0.748538,0.179867,1.169618,0.891611,0.937972,1.913545,1.265633,1.496658,0.347968,1.370786,1.623587,1.441953,0.349022,0.044401,0.433147,1.512506,0.908221,1.260697,1.179163,0.788047,2.634763,3.203953,1.413651,2.855107,2.117271,1.361689,0.682980,1.084290,2.651419,1.582776,2.098994,1.031491,0.974772,1.707198,1.686622,1.239945,0.910838,1.352740,0.491108,1.342217,1.562023,1.200723,0.186461,1.286620,0.674143,0.373987,1.223948,0.337008,0.731928,-0.281357,1.175176,0.739647,0.678546,1.123452,1.113558,0.726267,1.512785,2.043131,0.448006,2.183572,1.939738,1.332821,1.511216,1.018529,0.711897,2.307974,1.003607,1.043469,0.970110,1.871644,2.915320,2.135818,2.253222,2.031102,2.518522,0.921337,1.458060,0.792877,2.320386,1.245085,1.306543,0.180450,0.968659,0.164574,1.183990,1.127263,1.228676,0.268963,0.300094,0.112850,0.891135,1.174422,0.573210,0.442882,0.911176,0.576186,0.968842,0.497431,0.082573,0.328009,-0.056750,0.964732,0.475686,-0.015669,0.770434,-0.127846,0.252855,0.060786,0.259900,0.518966,0.403355,-0.108861,0.481694,0.843791,0.464014,0.797863,0.415165,-0.252576,0.253181,0.015900,0.816841,0.834992,1.134047,0.432330,-0.083651,0.403472,0.308811,0.595590,0.620263,0.476225,0.213429,0.151142,0.977679,0.609814,1.004169,0.006015,1.024273,-0.661786,0.966018,0.854617,1.302198,1.527325,0.915777,1.430393,0.901922,0.386083,1.077657,1.109994,1.011158,0.950192,0.705919,0.214606,0.504625,0.684422,2.340174,1.532737,0.335218,0.831516,0.047987,0.790953,0.180184,0.324308,-0.272377,0.905840,0.459948,0.725185,0.757486,0.601972,0.792880,0.643391,0.734002,0.462234,0.766695,0.584122,0.899299,-0.014713,0.345940,0.288591,0.433621,0.370425,-0.004842,0.715316,0.537479,0.353948,-0.070736,0.733348,2.617352,1.361348,1.910328,0.835367,0.582507,1.705017,0.016519,1.130763,1.239875,0.619161,1.271028,0.375051,1.618386,2.060941,2.005157,0.408412,1.442482,1.966379,1.826671,1.906063,1.787012,1.913379,1.512496,1.719286,1.540749,1.192583,0.891827,1.072452,1.310453,1.257718,0.361218,1.002642,0.679120,0.858933,0.483172,2.331961,0.769973,2.426531,1.043434,-0.013113,1.534898,0.414657,2.561237,1.034985,1.472958,1.820603,1.793376,0.750470,0.275797,0.158138,1.543065,3.001387,0.910785,0.312544,1.413810,0.682202,0.431473,0.699683,0.236938,0.660202,0.558848,0.770453,0.774522,0.302053,-0.014613,0.670073,0.680634,0.521799,0.413918,0.118622,0.054225,1.476402,0.750599,1.070832,0.656672,0.735283,0.655755,0.521329,1.149085,0.609177,1.329232,0.787796,0.590753,0.680317,0.482503,0.278492,0.274075,0.802026,-0.231484,0.459489,0.211343,0.638176,1.222268,1.115528,0.459551,1.168715,0.865882,2.308868,1.768827,1.772107,1.436916,2.629544,1.882291,1.167575,1.276729,1.290020,1.995226,1.019229,0.594557,1.522644,0.117587,1.206995,2.319513,1.735668,2.944532,1.240982,1.219626,0.480956,1.072602,1.563292,1.527690,1.554134,-0.171939,2.048136,2.266864,0.480223,0.377766,0.535421,0.483942,0.196052,0.054045,0.258394,0.798785,0.649261,0.534345,0.055019,0.243180,0.657302,0.706254,0.278796,0.615366,0.380890,0.276906,0.243272,0.551661,0.863483,0.440630,0.798439,0.586324,1.195738,1.009371,0.579874,0.372017,0.553026,1.332955,0.609213,0.748541,0.866370,0.833495,0.095667,0.401603,0.856719,0.345496,0.635498,-0.285572,0.581866,1.206264,0.733067,0.724146,0.666546,0.385326,0.415055,0.211203,0.151550,0.555357,0.802586,0.660676,0.582425,0.439438,0.007551,0.473846,0.635077,0.664830,0.365472,1.127987,0.490207,0.652832,0.161658,0.065514,1.235893,0.290362,0.615830,0.416137,0.363767,0.219042,0.478298,0.331287,0.944739,2.051951,1.434720,1.836504,1.286601,0.294812,0.638068,-0.268254,1.295887,0.480747,0.208821,0.855976,0.269729,0.914762,1.275772,1.534048,0.170300,0.904900,1.051581,2.511212,1.917051,2.286645,1.777414,0.811605,0.860505,0.603685,1.335006,1.213162,0.271023,0.757828,2.094239,2.078206,1.625170,0.176490,1.065432,1.349759,0.589513,1.268969,0.529612,0.806124,0.787966,0.234864,0.830407,0.542952,0.836920,0.882701,0.335129,0.472933,0.341546,0.973407,1.307088,0.238449,0.693853,0.270519,1.009405,1.406329,-0.570462,0.891382,2.428313,1.094265,1.783390,1.040052,0.861644,0.203140,0.636084,1.450647,0.847084,0.969826,1.414097,0.198779,0.687490,0.522458,0.053903,0.408947,0.569543,0.568064,1.260199,1.261880,0.220552,0.911870,0.127079,0.707255,0.832581,0.449130,0.520675,0.211952,0.028705,0.457395,0.659272,-0.242477,0.516395,0.274078,0.276373,-0.002293,1.556410,0.580424,1.316127,0.525453,1.171136,1.053243,0.751699,1.414437,1.256097,0.827122,0.528781,0.961972,1.363223,1.003883,1.148504,0.225404,1.195053,0.768401,0.313329,1.179661,0.773226,1.113219,0.382722,1.048375,0.141665,-0.123342,0.130314,0.641644,0.611059,0.549747,0.975831,0.504473,1.235451,0.285936,0.730503,0.538290,0.257493,1.136582,1.346161,0.541544,0.437374,1.060669,0.729838,0.636636,0.318465,0.451532,0.506974,0.615706,0.500159,0.489210,1.483199,1.510834,1.596820,1.491310,1.879318,2.831411,2.084957,1.944947,0.376177,0.515301,0.460145,0.879215,0.721696,1.080739,1.765914,1.471910,1.428790,1.086297,1.964416,0.615656,1.898981,-0.075818,1.583751,2.361304,1.504454,2.179671,1.393165,1.349864,0.577481,0.391314,0.677646,0.845260,0.228248,1.348949,2.012490,1.078521,0.864721,1.179847,1.032699,1.126812,0.288954,0.454042,1.477512,1.082741,1.608072,0.592619,2.113638,1.145185,1.615081,2.350135,1.685549,2.680238,1.630852,1.538801,1.061139,0.531765,0.705720,0.588122,0.045722,0.383626,0.244271,0.485260,-0.078997,0.236802,0.547583,0.478837,0.163961,0.260605,0.503110,-0.145127,0.175072,0.736469,0.674464,0.658005,0.052520,1.203791,0.542004,-0.346794,0.502682,-0.070346,-0.137591,0.480990,0.242443,0.538687,-0.016141,1.853662,0.738113,0.263812,1.025296,0.339119,2.151932,0.823091,1.874306,1.295417,0.359966,1.170363,1.245068,1.097720,0.956104,1.567615,1.988063,0.468615,1.263019,1.009524,1.185564,-0.006185,1.580293,0.845020,0.277008,0.424592,0.637632,1.007791,0.395600,0.732423,0.883847,2.054144,0.851393,0.624371,1.228377,0.652804,0.837658,0.812056,1.222096,0.358725,0.454917,0.649817,1.614692,0.465409,0.997623,3.101304,5.273570,3.414006,1.914391,1.120832,0.586741,0.785430,1.303540,0.849093,0.836327,1.705207,0.546375,0.711407,1.403707,0.856065,0.252782,0.369002,0.503240,1.007064,1.623281,0.924029,0.295286,0.406171,0.331450,1.161005,0.836686,1.225883,0.489773,0.859426,0.606825,1.117295,0.651693,0.679639,0.712851,0.868965,0.540321,0.166122,0.276101,0.001476,0.559463,0.407958,0.811634,0.092167,1.301193,-0.262239,0.949893,0.583211,0.584087,1.032614,1.276221,0.273619,4.095887,4.056882,3.837631,-0.006470,-0.087727,1.890815,1.914572,0.330307,0.973817,0.147949,0.041069,-0.033778,0.791540,1.289544,0.638887,1.669934,0.590356,0.676710,0.923502,0.965393,0.468341,0.393931,1.445475,0.640377,1.262432,0.965532,0.881629,0.369243,1.250458,0.873058,2.320135,1.536949,0.853499,2.091005,0.379653,3.375722,3.502891,0.384292,4.117616,0.458369,0.854213,0.486973,0.294892,0.679739,0.512834,0.304734,0.858806,1.422484,0.923356,1.084787,1.118456,1.256853,0.523900,0.891699,1.622718,1.991983,1.052135,1.113966,0.281608,0.257447,0.618557,-0.041693,0.146650,1.262610,1.184756,0.834426,0.412848,0.624662,0.293311,1.027782,0.998640,0.735239,0.645252,0.534118,0.528553,0.753163,0.362140,-0.153894,1.114657,0.725707,0.611680,0.727001,0.710794,0.343803,0.926777,0.814910,0.744707,3.921001,3.771756,3.615430,0.406582,0.267640,0.327475,1.804561,1.266996,0.074195,0.505212,1.241513,0.782211,0.621392,1.240858,1.032221,0.052604,0.162458,0.712631,0.721223,1.214213,1.313806,0.108857,0.318761,0.753520,-0.007022,0.131622,0.121098,0.448717,0.583963,0.462094,0.758642,0.476337,0.273397,0.615087,0.783056,0.320645,1.926426,2.567057,2.451532,0.125939,0.385409,-0.182078,0.936257,0.788181,0.284974,0.248804,0.197322,0.373591,1.422602,1.722105,1.889732,-0.031607,0.110257,0.392775,0.379231,0.397524,1.427561,0.483742,0.100910,-0.115655,1.921954,1.531023,1.765195,0.362692,0.642126,0.996490,1.351157,0.510810,0.583058,0.667169,-0.032673,1.148892,1.278514,0.100629,0.318408,2.162623,2.201889,2.460349,-0.024192,-0.064966,0.388374
    
])

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
ax.pcolormesh(
    np.rad2deg(ψ_space), 
    np.rad2deg(θ_space), 
    E_space.T
)
ax.set_xlabel("Azimuth ($\degree$)")
ax.set_ylabel("Elevation ($\degree$)")
plt.title("The Energy-Map from the Physical Prototype with a Metronome")
plt.savefig("./main_method/physical/map.png", dpi = 1024)
plt.show()
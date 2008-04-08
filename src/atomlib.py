"""
Form factors from Int. Tab. Cryst Sect. C 6.1.1.4

"""

formfactor = {\
    'H'  : [ 0.493000, 0.322910, 0.140190, 0.040810, 10.510910, 26.125730, 3.142360, 57.799770, .003038],
    'HE' : [ 0.873400, 0.630900, 0.311200, 0.178000, 9.103710, 3.356800, 22.927630, 0.982100, .00640],
    'LI' : [ 1.128200, 0.750800, 0.617500, 0.465300, 3.954600, 1.052400, 85.390580, 168.261200, .03770],
    'BE' : [ 1.591900, 1.127800, 0.539100, 0.702900, 43.642750, 1.862300, 103.483100, 0.542000, .03850],
    'B'  : [ 2.054500, 1.332600, 1.097900, 0.706800, 23.218520, 1.021000, 60.349870, 0.140300, 0.19320],
    'C'  : [ 2.310000, 1.020000, 1.588600, 0.865000, 20.843920, 10.207510, 0.568700, 51.651250, .21560],
    'N'  : [12.212610, 3.132200, 2.012500, 1.166300, 0.005700, 9.893310, 28.997540, 0.582600, 11.52901],
    'O'  : [ 3.048500, 2.286800, 1.546300, 0.867000, 13.277110, 5.701110, 0.323900, 32.908940, .25080],
    'F'  : [ 3.539200, 2.641200, 1.517000, 1.024300, 10.282510, 4.294400, 0.261500, 26.147630, .27760],
    'NE' : [ 3.955300, 3.112500, 1.454600, 1.125100, 8.404210, 3.426200, 0.230600, 21.718410, .35150],
    'NA' : [ 4.762600, 3.173600, 1.267400, 1.112800, 3.285000, 8.842210, 0.313600, 129.424100, .67600],
    'MG' : [ 5.420410, 2.173500, 1.226900, 2.307300, 2.827500, 79.261180, 0.380800, 7.193710, .85840],
    'AL' : [ 6.420210, 1.900200, 1.593600, 1.964600, 3.038700, 0.742600, 31.547240, 85.088680, .11510],
    'SI' : [ 6.291510, 3.035300, 1.989100, 1.541000, 2.438600, 32.333740, 0.678500, 81.693790, .14070],
    'P'  : [ 6.434510, 4.179100, 1.780000, 1.490800, 1.906700, 27.157040, 0.526000, 68.164570, .11490],
    'S'  : [ 6.905310, 5.203410, 1.437900, 1.586300, 1.467900, 22.215120, 0.253600, 56.172070, .86690],
    'CL' : [ 11.460410, 7.196410, 6.255610, 1.645500, 0.010400, 1.166200, 18.519420, 47.778460, 9.55741],
    'AR' : [ 7.484510, 6.772310, 0.653900, 1.644200, 0.907200, 14.840710, 43.898350, 33.392930, .44450],
    'K'  : [ 8.218610, 7.439810, 1.051900, 0.865900, 12.794910, 0.774800, 213.187200, 41.684160, .42280],
    'CA' : [ 8.626610, 7.387310, 1.589900, 1.021100, 10.442110, 0.659900, 85.748490, 178.437200, .37510],
    'SC' : [ 9.189010, 7.367910, 1.640900, 1.468000, 9.021310, 0.572900, 136.108100, 51.353150, .33290],
    'TI' : [ 9.759510, 7.355810, 1.699100, 1.902100, 7.850810, 0.500000, 35.633830, 116.105100, .28070],
    'V'  : [ 10.297110, 7.351110, 2.070300, 2.057100, 6.865710, 0.438500, 26.893830, 102.478100, .21990],
    'CR' : [ 10.640610, 7.353710, 3.324000, 1.492200, 6.103810, 0.392000, 20.262620, 98.739990, .18320],
    'MN' : [ 11.281910, 7.357310, 3.019300, 2.244100, 5.340910, 0.343200, 17.867420, 83.754380, .08960],
    'FE' : [ 11.769510, 7.357310, 3.522200, 2.304500, 4.761110, 0.307200, 15.353510, 76.880580, .03690],
    'CO' : [ 12.284110, 7.340910, 4.003400, 2.348800, 4.279100, 0.278400, 13.535910, 71.169270, .01180],
    'NI' : [ 12.837610, 7.292010, 4.443800, 2.380000, 3.878500, 0.256500, 12.176310, 66.342160, .03410],
    'CU' : [ 13.338010, 7.167610, 5.615810, 1.673500, 3.582800, 0.247000, 11.396610, 64.812670, .19100],
    'ZN' : [ 14.074310, 7.031810, 5.165210, 2.410000, 3.265500, 0.233300, 10.316310, 58.709760, .30410],
    'GA' : [ 15.235410, 6.700610, 4.359100, 2.962300, 3.066900, 0.241200, 10.780510, 61.413570, .71890],
    'GE' : [ 16.081620, 6.374710, 3.706800, 3.683000, 2.850900, 0.251600, 11.446810, 54.762560, .13130],
    'AS' : [ 16.672320, 6.070110, 3.431300, 4.277900, 2.634500, 0.264700, 12.947910, 47.797260, .53100],
    'SE' : [ 17.000630, 5.819610, 3.973100, 4.354300, 2.409800, 0.272600, 15.237210, 43.816350, .84090],
    'BR' : [ 17.178920, 5.235810, 5.637710, 3.985100, 2.172300, 16.579620, 0.260900, 41.432850, .95570],
    'KR' : [ 17.355510, 6.728610, 5.549310, 3.537500, 1.938400, 16.562320, 0.226100, 39.397230, .82500],
    'RB' : [ 17.178420, 9.643510, 5.139900, 1.529200, 1.788800, 17.315120, 0.274800, 164.934200, .48730],
    'SR' : [ 17.566310, 9.818410, 5.422000, 2.669400, 1.556400, 14.098810, 0.166400, 132.376100, .50640],
    'Y'  : [ 17.776020, 10.294610, 5.726300, 3.265880, 1.402900, 12.800610, 0.125600, 104.354100, .91213],
    'ZR' : [ 17.876530, 10.948010, 5.417330, 3.657210, 1.276180, 11.916010, 0.117620, 87.662780, .06929],
    'NB' : [ 17.614230, 12.014410, 4.041830, 3.533460, 1.188650, 11.766010, 0.204790, 69.795760, .75591],
    'MO' : [ 3.702500, 17.235630, 12.887610, 3.742900, 0.277200, 1.095800, 11.004010, 61.658460, .38750],
    'TC' : [ 19.130130, 11.094810, 4.649020, 2.712630, 0.864130, 8.144880, 21.570720, 86.847270, .40429],
    'RU' : [ 19.267430, 12.918210, 4.863370, 1.567560, 0.808520, 8.434680, 24.799740, 94.292890, .37875],
    'RH' : [ 19.295720, 14.350110, 4.734250, 1.289180, 0.751540, 8.217590, 25.874940, 98.606290, .32800],
    'PD' : [ 19.331920, 15.501720, 5.295370, 0.605840, 0.698660, 7.989300, 25.205230, 76.898680, .26593],
    'AG' : [ 19.280820, 16.688520, 4.804510, 1.046300, 0.644600, 7.472610, 24.660540, 99.815700, .17900],
    'CD' : [ 19.221420, 17.644420, 4.461000, 1.602900, 0.594600, 6.908910, 24.700840, 87.482570, .06941],
    'IN' : [ 19.162410, 18.559620, 4.294800, 2.039600, 0.547600, 6.377610, 25.849930, 92.802990, .93911],
    'SN' : [ 19.188920, 19.100520, 4.458500, 2.466300, 5.830310, 0.503100, 26.890930, 83.957180, .78211],
    'SB' : [ 19.641820, 19.045520, 5.037110, 2.682700, 5.303400, 0.460700, 27.907440, 75.282580, .59091],
    'TE' : [ 19.964420, 19.013820, 6.144880, 2.523900, 4.817420, 0.420890, 28.528440, 70.840360, .35200],
    'I'  : [ 20.147220, 18.994920, 7.513810, 2.273500, 4.347000, 0.381400, 27.766040, 66.877670, .07120],
    'XE' : [ 20.293320, 19.029820, 8.976710, 1.990000, 3.928200, 0.344000, 26.465940, 64.265870, .71180],
    'CS' : [ 20.389220, 19.106220, 10.662010, 1.495300, 3.569000, 0.310700, 24.387940, 213.904200, .33520],
    'BA' : [ 20.336120, 19.297030, 10.888010, 2.695900, 3.216000, 0.275600, 20.207320, 167.202200, .77310],
    'LA' : [ 20.578020, 19.599010, 11.372710, 3.287190, 2.948170, 0.244480, 18.772610, 133.124100, .14678],
    'CE' : [ 21.167110, 19.769520, 11.851310, 3.330490, 2.812190, 0.226840, 17.608320, 127.113100, .86264],
    'PR' : [ 22.044020, 19.669720, 12.385610, 2.824280, 2.773930, 0.222090, 16.766920, 143.644100, .05830],
    'ND' : [ 22.684520, 19.684720, 12.774010, 2.851370, 2.662480, 0.210630, 15.885020, 137.903100, .98486],
    'PM' : [ 23.340520, 19.609530, 13.123510, 2.875160, 2.562700, 0.202090, 15.100910, 132.721100, .02876],
    'SM' : [ 24.004240, 19.425830, 13.439610, 2.896040, 2.472740, 0.196450, 14.399610, 128.007100, .20963],
    'EU' : [ 24.627440, 19.088620, 13.760310, 2.922700, 2.387900, 0.194200, 13.754610, 123.174100, .57450],
    'GD' : [ 25.070940, 19.079820, 13.851810, 3.545450, 2.253410, 0.181950, 12.933110, 101.398100, .41960],
    'TB' : [ 25.897630, 18.218520, 14.316710, 2.953540, 2.242560, 0.196140, 12.664810, 115.362100, .58324],
    'DY' : [ 26.507030, 17.638320, 14.559620, 2.965770, 2.180200, 0.202170, 12.189910, 111.874100, .29728],
    'HO' : [ 26.904940, 17.294020, 14.558310, 3.638370, 2.070510, 0.197940, 11.440710, 92.656690, .56797],
    'ER' : [ 27.656340, 16.428530, 14.977910, 2.982330, 2.073560, 0.223550, 11.360410, 105.703100, .92047],
    'TM' : [ 28.181930, 15.885120, 15.154210, 2.987060, 2.028590, 0.238850, 10.997510, 102.961100, .75622],
    'YB' : [ 28.664140, 15.434510, 15.308710, 2.989630, 1.988900, 0.257120, 10.664710, 100.417100, .56673],
    'LU' : [ 28.947630, 15.220810, 15.100010, 3.716010, 1.901820, 9.985200, 0.261030, 84.329880, .97629],
    'HF' : [ 29.144040, 15.172610, 14.758610, 4.300130, 1.832620, 9.599910, 0.275120, 72.029080, .58155],
    'TA' : [ 29.202440, 15.229310, 14.513510, 4.764920, 1.773330, 9.370470, 0.295980, 63.364470, .24355],
    'W'  : [ 29.081830, 15.430010, 14.432710, 5.119830, 1.720290, 9.225910, 0.321700, 57.056060, .88751],
    'RE' : [ 28.762130, 15.718920, 14.556410, 5.441740, 1.671910, 9.092280, 0.350500, 52.086150, 0.47201],
    'OS' : [ 28.189440, 16.155010, 14.930510, 5.675900, 1.629030, 8.979490, 0.382660, 48.164750, 1.00051],
    'IR' : [ 27.304930, 16.729610, 15.611520, 5.833780, 1.592790, 8.865540, 0.417920, 45.001140, 1.47221],
    'PT' : [ 27.005940, 17.763920, 15.713120, 5.783710, 1.512930, 8.811750, 0.424590, 38.610340, 1.68831],
    'AU' : [ 16.881930, 18.591320, 25.558240, 5.860010, 0.461100, 8.621610, 1.482600, 36.395630, 2.06581],
    'HG' : [ 20.680920, 19.041720, 21.657520, 5.967610, 0.545000, 8.448410, 1.572900, 38.324630, 2.60891],
    'TL' : [ 27.544630, 19.158420, 15.538020, 5.525940, 0.655150, 8.707520, 1.963470, 45.814960, 3.17461],
    'PB' : [ 31.061740, 13.063710, 18.442020, 5.969610, 0.690200, 2.357600, 8.618010, 47.257950, 3.41181],
    'BI' : [ 33.368940, 12.951010, 16.587720, 6.469210, 0.704000, 2.923800, 8.793710, 48.009350, 3.57821],
    'PO' : [ 34.672640, 15.473310, 13.113810, 7.025890, 0.701000, 3.550780, 9.556430, 47.004550, 3.67701],
    'AT' : [ 35.316330, 19.021120, 9.498880, 7.425190, 0.685870, 3.974580, 11.382410, 45.471560, 3.71081],
    'RN' : [ 35.563140, 21.281620, 8.003710, 7.443310, 0.663100, 4.069100, 14.042210, 44.247340, 3.69051],
    'FR' : [ 35.929930, 23.054720, 12.143910, 2.112530, 0.646450, 4.176190, 23.105220, 150.645100, 3.72471],
    'RA' : [ 35.763030, 22.906420, 12.473910, 3.210970, 0.616340, 3.871350, 19.988720, 142.325100, 3.62111],
    'AC' : [ 35.659730, 23.103230, 12.597710, 4.086550, 0.589090, 3.651550, 18.599010, 117.020100, 3.52661],
    'TH' : [ 35.564530, 23.421920, 12.747310, 4.807040, 0.563360, 3.462040, 17.830920, 99.172300, 3.43141],
    'PA' : [ 35.884740, 23.294820, 14.189110, 4.172870, 0.547750, 3.415190, 16.923520, 105.251100, 3.42871],
    'U'  : [ 36.022840, 23.412830, 14.949110, 4.188000, 0.529300, 3.325300, 16.092730, 100.613100, 3.39661],
    'NP' : [ 36.187440, 23.596420, 15.640220, 4.185500, 0.511930, 3.253960, 15.362220, 97.490890, 3.35731],
    'PU' : [ 36.525440, 23.808320, 16.770720, 3.479470, 0.499380, 3.263710, 14.945510, 105.980100, 3.38121]}


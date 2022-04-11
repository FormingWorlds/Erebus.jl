using BenchmarkTools
using TimerOutputs

using HydrologyPlanetesimals

const to = TimerOutput()


# initialize parameters
const static_parameters = StaticParameters(
    Nx = 141,
    Ny = 141,
    Nxmc = 4,
    Nymc = 4
)
# initialize dynamic parameters from static parameters
const dynamic_parameters = DynamicParameters(static_parameters)
# const basicnodes = BasicNodes(static_parameters)
# const vxnodes = VxNodes(static_parameters)
# const vynodes = VyNodes(static_parameters)
# const pnodes = PNodes(static_parameters)
const markers = MarkerArrays(static_parameters.startmarknum)
initmarkers!(markers, static_parameters, dynamic_parameters)
import numpy as np
import math
import random

class Weed(object):

    # Attributes
    minSeeds = 0
    maxSeeds = 2

    xInfo = dict()




#
# properties
# x % induction
# fitness % objective
# function
# value
# objFunc
# end
#
# properties(SetAccess=private)
# minSeeds = 0
# maxSeeds = 2
#
# xInfo = struct();
# xDim
# x_min
# x_max
# end
#
# methods

    def __init__(self, x, xInfo):
        self.x = x  # decision variable

        self.xInfo = xInfo
        self.xDim = xInfo['xDim']
        self.x_min = xInfo['x_min']
        self.x_max = xInfo['x_max']

        self.evaluateFitness()

# function
# obj = Weed(x, xInfo, objFunc)
#
# if nargin > 0
#     obj.x = x;
#
#     obj.xInfo = xInfo;
#     obj.xDim = xInfo.xDim;
#     obj.x_min = xInfo.x_min;
#     obj.x_max = xInfo.x_max;
#
#     obj.objFunc = objFunc;
#
#     obj.evaluateFitness()
# end
# end
#

    def evaluateFitness(self):
        self.fitness = 1            # define objective function here

# function
# evaluateFitness(obj)
# obj.fitness = obj.objFunc(obj.x);
# end
#

    def createSeeds(self, minFit, maxFit, sig_curr):
        n_seeds = int(math.floor(self.minSeeds + (self.fitness - minFit)* (self.maxSeeds - self.minSeeds) / (maxFit - minFit)))

        seeds = []

        if n_seeds > 0:
            for i in range(n_seeds):
                x_new = self.spatialDisp(sig_curr)
                seeds.append(Weed(x_new, self.xInfo))


# function[seeds, n_seeds] = createSeeds(obj, minFit, maxFit, sig_curr)
# n_seeds = floor(obj.minSeeds + (obj.fitness - minFit)...
#                 * (obj.maxSeeds - obj.minSeeds) / (maxFit - minFit));
#
# if n_seeds == 0
#     seeds = [];
# else
#
#     seeds(n_seeds, 1) = Weed();
#
#     for i=1:n_seeds
#     x_new = spatialDisp(obj, sig_curr);
#     seeds(i, 1) = Weed(x_new, obj.xInfo, obj.objFunc);
# end
# end
# end
#

    def spatialDisp(self, sig_curr):
        x_new = np.zeros(self.xDim)

        for i in np.arange(0,self.xDim,1):
            xi_new = self.x[i] + random.gauss(0, sig_curr)

            if xi_new > self.x_max[i]:
                xi_new = self.x_max[i]

            if xi_new < self.x_min[i]:
                xi_new = self.x_min[i]

            x_new[i] = xi_new

# function
# x_new = spatialDisp(obj, sig_curr)
# x_new = zeros(obj.xDim, 1);
#
# for i=1:obj.xDim
# xi_new = obj.x(i) + sig_curr * randn;
#
# if xi_new > obj.x_max(i)
#     xi_new = obj.x_max(i);
# end
#
# if xi_new < obj.x_min(i)
#     xi_new = obj.x_min(i);
# end
#
# x_new(i) = xi_new;
# end
#
# end
#
# end
#
# end

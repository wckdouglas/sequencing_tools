import numpy as np
from ..utils import SeqUtilsError
import logging
logging.basicConfig(level = logging.INFO)
logger = logging.getLogger('Regression')


class Bootstrap:
    def __init__(self, seed=123):
        '''
        boostrap 1d array
        usage:
        xs = np.arange(100)
        bs = Bootstrap(seed=123)
        for idx in bs.bootstrap(xs, group_size=50, n_boots=10):
            print(xs[idx].mean())
        '''
        self.rng = np.random.RandomState(seed)

    def bootstrap(self, xs, group_size=100, n_boots = 100):
        '''
        input:
            xs: 1d np.array
            group_size: number of values in each bootstrap iteration
            n_boots: how many bootstrap groups
        output:
            iterator: bootstrapped
        '''
        xs = np.array(xs)
        total_size = xs.shape[0]
        logger.info('Total size for bootstrap: %i' %total_size)
        if group_size > total_size:
            raise SeqUtilsError('Group size > input array size')
    
        for i in range(n_boots):
            idx = self.rng.randint(0, total_size, group_size)
            yield idx


class GradientDescent():
    def __init__(self, lr = 0.001, max_iter = 10000, limit = 1e-4, verbose = False):
        '''
        A stochastic gradient descent model using Adam optimizer 
        for solving linear regression with no intercept:
        y = bx . solving for b
        input:
            lr: learning rate for how gradient should change the weight in each iteration
            max_iter: how many iteration to run if not converging
            limit: how small of a difference between iteration we would tolerate as converge?
        test:
        import logging
        import numpy as np
        logging.basicConfig(level=logging.INFO)    
        np.random.seed(123)
        X = 2 * np.random.rand(100)
        y = 3 * X+np.random.randn(100)
        gd = GradientDescent()
        gd.fit(X,y)
        plt.plot(gd.losses)
        '''

        # static parameter 
        self.learning_rate = lr
        self.max_iter = max_iter
        self.verbose = verbose
        self.print = self.max_iter // 5
        self.b = np.random.rand()
        self.diff = 1
        self.beta_1 = 0.9
        self.beta_2 = 0.999
        self.initialize()
        self.epsilon = 1e-8 # avoid division by zero
        self.limit = limit
        self.bootstrap = Bootstrap(seed=123)
        self.logger = logging.getLogger('Gradient Descent')
        self.initialize()

    def initialize(self):
        # initialize parameters
        self.moving_average_gradient = 0
        self.moving_average_squared_gradient = 0
        self._iter = 0
        self.x = None
        self.y = None
        self.n = 0
        self.cost = 0
        self.losses = np.zeros(int(self.max_iter))
        self.gradients = np.zeros(int(self.max_iter))
        self.diffs = np.zeros(int(self.max_iter))
        if self.verbose:
            self.logger.info('Initialized parameters')
        self.converge = False

    
    def Cost(self):
        '''
        compute error with new regression coefficien
        '''
        self.cost = self.y - self.x * self.b  # error = y - predicted_y 
        self.gradient = np.median( self.x * self.cost ) # using a median cost
        self.losses[self._iter - 1] = np.abs(self.gradient)
        self.gradients[self._iter - 1] = self.gradient
    
    
    def Adam_update(self):
        '''
        Adam optimizer
        https://github.com/sagarvegad/Adam-optimizer/blob/master/Adam.py
        '''
        self.moving_average_gradient = self.beta_1*self.moving_average_gradient + (1-self.beta_1)*self.gradient
        self.moving_average_squared_gradient = self.beta_2*self.moving_average_squared_gradient + (1-self.beta_2)*(self.gradient * self.gradient)
        m_cap = self.moving_average_gradient / (1-(self.beta_1**self._iter)) #calculates the bias-corrected estimates
        v_cap = self.moving_average_squared_gradient / (1-(self.beta_2**self._iter)) #calculates the bias-corrected estimates

        self.diff = self.learning_rate * m_cap / (np.sqrt(v_cap) + self.epsilon)
        self.diffs[self._iter - 1] = self.diff
        self.b -=  self.diff
        self.diff = np.abs(self.diff)

    def fit(self, x, y):
        '''
        fitting B for 
        
        y = Bx
        '''
        self._iter = 0
        self.n = len(x)
        bootstrap_idx = self.bootstrap.bootstrap(x, group_size=len(x)//5, 
                                                n_boots=int(self.max_iter))

        assert(x.ndim == 1)
        assert(x.ndim == y.ndim)

        while self._iter < self.max_iter and self.diff < self.epsilon:
            idx = next(bootstrap_idx)
            self.x = x[idx]
            self.y = y[idx]
            self._iter += 1
            self.Cost()
            self.Adam_update()
            if (self._iter % self.print == 0 or self._iter == 1) and self.verbose:
                self.logger.info('%i iteration: gradient %.3f; b %.7f; diff %.7f' %(self._iter, self.gradient, self.b, self.diff))
        
        if self.diff > self.limit:
            if self.verbose:
                self.logger.warning('b is not converged, please consider increasing max_iter')
        elif self.verbose:
            self.converge = True
            self.logger.info('Converged at the %ith iteration: gradient %.3f; b %.7f; diff %.7f' %(self._iter, self.gradient, self.b, self.diff))
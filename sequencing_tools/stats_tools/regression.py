import numpy as np
import logging
from ..utils import SeqUtilsError
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
            #raise SeqUtilsError('Group size > input array size')
            raise ValueError('Group size > input array size')
    
        for i in range(n_boots):
            idx = self.rng.randint(0, total_size, group_size)
            yield idx
class GradientDescent():
    def __init__(self, 
                 lr = 0.001, max_iter = 10000, 
                 limit = 1e-4, verbose = False,
                 seed = 123,
                method = 'mean'):
        '''
        A stochastic gradient descent model using Adam optimizer 
        for solving linear regression with no intercept:
        y = bx . solving for b
        input:
            lr: learning rate for how gradient should change the weight in each iteration
            max_iter: how many iteration to run if not converging
            limit: how small of a difference between iteration we would tolerate as converge?
            method: mean or median
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
        assert(method in ['mean', 'median'])
        self.method = method
        self.learning_rate = lr
        self.max_iter = int(max_iter)
        self.verbose = verbose
        self.print = self.max_iter // 5
        self.B = None # regression coefficients
        self.diff = 1
        self.beta_1 = 0.9
        self.beta_2 = 0.999
        self.epsilon = 1e-8 # avoid division by zero
        self.limit = limit
        self.bootstrap = Bootstrap(seed=seed)
        self.rng = np.random.RandomState(seed)
        self.logger = logging.getLogger('Gradient Descent')
        self._initialize()

    def _initialize(self):
        # initialize parameters
        self.moving_average_gradient = 0
        self.moving_average_squared_gradient = 0
        self._iter = 0
        self.X = None
        self.y = None
        self.n = 0
        self.n_coefficients = 0
        self.B = []
        self.B_history = []
        self.cost = 0
        self.last_cost = 0
        self.losses = []
        self.gradients = []
        self.gradient = []
        self.diffs = []
        self.diff = 0 
        if self.verbose:
            self.logger.info('Initialized parameters')
        self.converge = False
        if self.method == 'mean':
            self.summary = np.mean
        if self.method == 'median':
            self.summary = np.median


    
    def Cost(self):
        '''
        compute error with new regression coefficien
        '''
        self.residuals = self.y - np.matmul(self.X, self.B)  # error = y - predicted_y 
        self.cost = np.sqrt(self.summary(self.residuals**2))
        self.losses[self._iter - 1] = self.cost

        for i in range(self.n_coefficients):
            self.gradient[i] = self.summary( - self.X[:,i] * self.residuals[i] )
            self.gradients[self._iter - 1, i] = self.gradient[i]

    
    def Adam_update(self):
        '''
        Adam optimizer
        https://github.com/sagarvegad/Adam-optimizer/blob/master/Adam.py
        '''

        for i in range(self.n_coefficients):
            self.moving_average_gradient = self.beta_1*self.moving_average_gradient + (1-self.beta_1)*self.gradient[i]
            self.moving_average_squared_gradient = self.beta_2*self.moving_average_squared_gradient + (1-self.beta_2)*(self.gradient[i] * self.gradient[i])
            m_cap = self.moving_average_gradient / (1-(self.beta_1**self._iter)) #calculates the bias-corrected estimates
            v_cap = self.moving_average_squared_gradient / (1-(self.beta_2**self._iter)) #calculates the bias-corrected estimates

            self.diff = self.learning_rate * m_cap / (np.sqrt(v_cap) + self.epsilon)
            self.diffs[i] = self.diff
            self.B[i] -=  self.diff
        self.diffs = np.abs(self.diffs)

    def fit(self, X, y):
        '''
        fitting B for 
        
        y = Bx
        '''
        if X.ndim != 2:
            #raise SeqUtilsError("X must be 2 dimentional: do X.reshape(-1,1) if it's 1-d")
            raise ValueError("X must be 2 dimentional: do X.reshape(-1,1) if it's 1-d")

        self.orig_X = X
        self.orig_y = y
        self._iter = 0
        self.n = len(X)
        self.n_coefficients = X.shape[1]
        self.B = self.rng.rand(self.n_coefficients)
        self.B_history = np.zeros((self.max_iter, self.n_coefficients))
        self.gradient = np.zeros(self.n_coefficients)
        self.diffs = np.zeros(self.n_coefficients) + 1000
        self.losses = np.zeros(self.max_iter)
        self.gradients = np.zeros((self.max_iter, self.n_coefficients))
        self.bootstrap_idx = self.bootstrap.bootstrap(X, group_size=len(X)//10, 
                                                n_boots=int(self.max_iter))

        self._fit()
        self.logger.info('%i iteration: Cost %.3f; Diff %.7f' %(self._iter, self.cost, self.diffs.max()))
        while self._iter < self.max_iter and self.diffs.max() > self.limit:
            self._fit()
            if self._iter % self.print == 0  and self.verbose:
                self.logger.info('%i iteration: Cost %.3f, Diff %.7f' %(self._iter, self.cost, self.diffs.max()))
        
        if self.diffs.max() > self.limit:
            if self.verbose:
                self.logger.warning('b is not converged, please consider increasing max_iter')
        elif self.verbose:
            self.converge = True
            self.logger.info('Converged at the %ith iteration: Cost %.3f, Diff %.7f' %(self._iter, self.cost, self.diffs.max()))
            
    def _fit(self):
        self.last_B = self.B
        idx = next(self.bootstrap_idx)
        self.X, self.y = self.orig_X[idx], self.orig_y[idx]
        self._iter += 1
        self.Cost()
        self.Adam_update()
        self.B_history[self._iter - 1] = self.B
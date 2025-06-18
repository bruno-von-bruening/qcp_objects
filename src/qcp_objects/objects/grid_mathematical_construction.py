""" This function is for the mathematical constrcuction of grids!"""
from . import *

class math_radial_grid_base(BaseModel):
    class Config:
        arbitrary_types_allowed=True
    radial_points   : mynpar = np.array([])
    weights         : mynpar = np.array([])
    number_of_points : int=200

    def get_grid(self):
        return self.radial_points, self.weights
    def integrate(self, function_values:mynpar):
        assert function_values.shape==(self.number_of_points,)
        return np.sum(self.weights*function_values, axis=0)
    def test(self):
        @val_call
        def my_exp(input:mynpar|None=None, alpha:float=1., analytical:bool=False):
            norm=(8*np.pi)*np.pow(1/alpha, 3)
            if not analytical:
                assert isinstance(input, np.ndarray)
                return np.exp(-alpha*input)/norm
            else:
                return 1
        @val_call
        def my_gau(input:mynpar|None=None, alpha:float=1., analytical:bool=False):
            norm=np.pow( np.pi, 3/2 )/np.pow( alpha, 3/2)
            if not analytical:
                assert isinstance(input, np.ndarray)
                return np.exp(-alpha*input**2) / norm
            else:
                return 1.

        test_functions={
            'exp':my_exp,
            'gau':my_gau,
        }
        tests=[
            ('exp',[float(np.pow(10.,fac))]) for fac in range(-5,5) 
        ]+[
            ('gau',[float(np.pow(10.,fac))]) for fac in range(-5,5)
        ]
        for t in tests:
            func=test_functions[t[0]]
            args=[] if len(t)==1 else t[1]

            function_values=func(self.radial_points, *args)
            numerical_integral=self.integrate(function_values)
            analytic_integral=func(self.radial_points, *args, analytical=True)

            diff=numerical_integral-analytic_integral
            thres=1.e-10
            test_failed=  np.abs(diff)>=thres
            if test_failed: print(t, diff)


    # General tests

class Muray_Handy_Lamind_Euler_Maclaurin_Grid(math_radial_grid_base):
    """ Transform a series"""
    contraction     : int   =2
    scaling_factor  : float =3.
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        radial_points, weights=self.radial_construction(self.number_of_points)
        self.radial_points  = radial_points
        self.weights         = weights
    
    def radial_construction(self, number_of_grid_points=200):
        # This should be element specific (at a point of low high ->)

        def radial_transformation(number_of_grid_points):
            # First point and last point can be omitted since it is zero anway (i->i+1)
            # Hence only number_of_grid_points-2 points

            def core(n):
                """ (n/N-n)^m 
                m: contraction coefficient
                N: Total number of points
                n: 0<= n < N"""
                n=float(n)
                N=float(number_of_grid_points)
                m=self.contraction
                return np.pow(  
                    (n)/(N-n)
                ,  m  )
            
            grid_points=[ float(n)+.5 for n in range(number_of_grid_points) ]
            return np.array([ core(n)  for n in grid_points])
        
        def the_weights(number_of_grid_points):
            """ 
            Assume the transformation r-> r(n) then an radial symmetric integral becomes:
                int_0^inf 4pi r^2 rho(r) dr =  4pi * int_0^N r(n)^2  [ d(rho(r(n))) / dn ] dn 
            See the definition of r(n) in your function and calculte the derivative of n then one obtains:
                int ... = int_0^N   rho(r(n)) * w(r(n)) dn
            with the weightsZ
                w(r(n)) = 4pi r(n)^2 d(r(n))/dn
                        = ( 4pi alf^(3m) ) * (n)^(3m-1) / (N-n)^(3*m+1)
            the Euler maclaurin formula then lets us evaluate the integral as sum of mid points
                int_0^n f(n)w(n) dn = sum_n=0,N-1 f(n+.5)w(n+.5) 
            with f(n+.5) = rho( r(n+.5)) and w(n)=w(r(n+.5)) in our case 
            """

            def core(n:float):
                m=self.contraction
                N=float(number_of_grid_points)
                return N*m*np.pow(n, 3*m -1) / np.pow(N-n, 3*m+1)
            grid_points=[ float(n)+.5 for n in range(number_of_grid_points) ]
            return np.array([core(n) for n in grid_points ])

        alpha=self.scaling_factor
        radial_points=alpha*radial_transformation(number_of_grid_points)
        radial_weights=4.*np.pi*alpha**3 * the_weights(number_of_grid_points)
        
        return radial_points, radial_weights

class math_angular_grid(BaseModel):
    spherical_points:  List[Tuple[float,float,float]] # Points on unit sphere that describe angular grid
    spherical_weights: List[float]|None=None
    @property
    def num_points(self):
        return len(self.spherical_points)
    def __init__(self, *args, **kwargs): # automatically make weights (they are implicit if spherical symmetry is followed)
        super().__init__(*args, **kwargs)
        self.spherical_weights=np.ones(self.num_points) * 1./float(self.num_points)
        assert np.abs( np.sum(self.spherical_weights)- 1. )<1.e-15
    def get_grid(self):
        return np.array(self.spherical_points), np.array(self.spherical_weights)
    

class AtomicGrid(BaseModel):
    """ given a radial grid, angular grid(s) and an option for merging them (could have different angular grids for radial dependence) """
    class Config:
        arbitrary_types_allowed=True
    radial_grid     : Muray_Handy_Lamind_Euler_Maclaurin_Grid # allow for different ones here
    angular_grids   : List[math_angular_grid]|math_angular_grid
    angular_pruning : None = None 
    # Multiple different grids may be used (drop that for now) 

    # Auto generate in init (could also define them as property)
    shell_indices  : List[ Tuple[float,float] ] =[] # delimiters of every shell (needed for spherical averaging) use as indices=range(start, stop)
    _base_grid_points     : List[Tuple[float,float,float]] =[]
    _base_grid_weights    : List[float]                    =[]
    _grid_points            : List[Tuple[float,float,float]]|None =None
    _grid_weights       : List[float]|None                    =None
    center              : Tuple[float,float,float]|None          =None

    @property
    def num_radial_points(self):
        return self.radial_grid.number_of_points
    @property
    def grid_weights(self):
        if self._grid_weights is None: 
            assert not self._base_grid_weights is None
            self._grid_weights=self._base_grid_weights
        return self._grid_weights
    @property
    def grid_points(self):
        if self.center is None:
            return self._base_grid_points
        else:
            return self._base_grid_points+np.array(self.center)
    

    def __init__(self,*args, **kwargs):
        super().__init__(*args, **kwargs)
        self._base_grid_points, self._base_grid_weights, self.shell_indices = self.combine_grids()

    def combine_grids(self):
        new_points=[]
        new_weights=[]
        shell_indices=[]
        stop= 0
        rad_pnt, rad_wgt = self.radial_grid.get_grid()
        for i in range(self.num_radial_points):

            # Make a function as soon as this becomes more complicated! (liost of coords and list of float)
            if self.angular_pruning is None:
                if isinstance(self.angular_grids, math_angular_grid):
                    anggri=self.angular_grids
                else:
                    anggri=self.angular_grids[0]
                the_shell, shell_wgts = anggri.get_grid()
            else: raise Exception(f"Did not implement angular pruning yet")
            # Get the radial part (floats)

            # combine them
            new_points      += (rad_pnt[i]*the_shell).tolist()
            new_weights     += (rad_wgt[i]*shell_wgts).tolist()

            start=stop
            stop+=len(shell_wgts)
            shell_indices   += [(start, stop)]
        assert np.sum( np.array(shell_indices)[:,1]-np.array(shell_indices)[:,0])==len(new_weights), shell_indices
        new_points, new_weights, shell_indices = [ np.array(x) for x in [new_points, new_weights, shell_indices]]
        return new_points, new_weights, shell_indices

    def spherical_average(self,array):
        assert array.shape==self.grid_weights.shape
        
        ar=[]
        for lower,upper in self.shell_indices:
            average=np.average(array[lower:upper])
            ar.append(average)
        return np.array(ar)
    @val_call
    def integrate(self,evaluation:np.ndarray):
        assert evaluation.shape==self.grid_weights.shape, f"{evaluation.shape} {self.grid_weights.shape}"
        return float(np.sum(evaluation*self.grid_weights, axis=0))
    @val_call
    def recenter(self, new_origin:Tuple[float,float,float]|None=None, translation:Tuple[float,float,float]|None=None):
        if new_origin is not None:
            assert translation is None, f"Do not provide new origin and transation at same time"
            self.center=new_origin
        elif translation is not None:
            self.center=np.array(self.center)+np.array(translation)
        else: raise Exception(f"Provide either new origin or translation key")
    @val_call
    def embed_radial(self, values: np.ndarray):
        assert values.shape == self.radial_grid.weights.shape
        new_values=[]
        for i,(start,stop) in enumerate(self.shell_indices):
            new_values+= (values[i] * np.ones(stop-start)).tolist()
        return np.array(new_values)
        



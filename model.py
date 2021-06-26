from pyomo.environ import Objective, Constraint, Var, Set, Param
from pyomo.environ import NonNegativeReals
from pyomo.environ import SolverFactory, ConcreteModel, BuildAction
from pyomo.opt import SolverStatus, TerminationCondition


class whatif():
    """ """
    def __init__(self):
        self.model = ConcreteModel()
        self.input_data = {}
        initialise(self.model, input_data)
        rules = whatif_rules()


def initialise(self, m, input_data, rules):
    """Initialise sets to what-if."""
    sets = input_data['sets']
    R = m.R = Set(dimen=1, initialize=sets['R'])                    # Regions (e.g. countries)
    C = m.C = Set(dimen=1, initialize=sets['C'])                    # Consumers (e.g. remote, rural, urban)
    G = m.G = Set(dimen=1, initialize=sets['G'])                    # Goods (e.g. elect, maize, wheat, sugar, ...)
    A = m.A = Set(dimen=1, initialize=sets['A'])                    # Assets (e.g. turbine1_100MW_rogun, irrig100ha, )
    T = m.T = Set(dimen=1, initialize=sets['T'])                    # Time tuples (e.g. (y2020,m01,day) )
    N = m.N = Set(dimen=1, initialize=sets['N'])                    # Nodes (e.g. upstream, midstream, downstream)
    M = m.M = Set(dimen=1, initialize=sets['X'])                    # Mix of crops (e.g. wheat_cotton, ...)

    subsets = input_data['subsets']
    TY = m.TY = Set(dimen=1, within=T, initialize=subsets['TY'])    # Time-years (e.g. y2020,mAll,hAll)
    TM = m.TM = Set(dimen=1, within=T, initialize=subsets['TM'])    # Time-months (e.g. y2020,m012,hAll)
    TH = m.TH = Set(dimen=1, within=T, initialize=subsets['TH'])    # Time-hours (e.g. y2020,m012,hDay)

    # Variables - market trade of goods
    m.S = Var(A,G,TY, within=NonNegativeReals)                       # Supply of goods from asset at time
    m.D = Var(R,C,G,TY, within=NonNegativeReals)                     # Demand regional consumers for goods at time
    m.I = Var(A,G,TY, within=NonNegativeReals)                       # Import of goods by asset and time
    m.X = Var(A,G,TY, within=NonNegativeReals)                       # Export of goods by asset and time
    # Variables - water flows and use
    m.F = Var(N,TM, within=NonNegativeReals)                         # Outflow from node at time
    m.V = Var(A,TM, within=NonNegativeReals)                         # Volume of reservoir by reservoir asset and time
    m.U = Var(A,TM, within=NonNegativeReals)                         # Use of water (gross) by asset and time
    m.R = Var(A,TM, within=NonNegativeReals)                         # Return of water (gross) by asset and time
    # Variables - agriculture
    m.L = Var(A,M,TY, within=NonNegativeReals)                       # Land use by asset, crop mix and time
    # Variables - capital formation
    m.K = Var(A,TY, within=NonNegativeReals)                         # Capital formation, investment in time-year


class whatif_rules():
    """Rules for constraints and objective"""

    def R_balance_goods(self,model,r,g,t):
        """Goods market balance rule."""
        m = model
        sup = sum(m.S[a,g,t] for a in m.A_R['S',r])                 # Supply of good from assets by time
        exp = sum(m.X[a,g,t] for a in m.A_R['X',r])                 # Export along trade of good asset by time
        imp = sum(m.I[a,g,t] for a in m.A_R['I',r])                 # Import along trade of good asset by time
        dem = sum(m.D[r,c,g,t] for c in m.C)                        # Demand by consumer for good by time
        return sup + imp == exp + dem

    def R_balance_water(self,model,n,t):
        """Node water balance."""
        m = model
        ups = sum(m.F[nu,t] for nu in m.N_N[n])                     # Flow from upstream by time
        use = sum(m.U[a,t] for a in m.A_N['U',n])                   # Use of water by asset from node by time
        ret = sum(m.R[a,t] for a in m.A_N['R',n])                   # Return of water by asset from node by time
        return m.F[n,t] == (ups + ret - use)*(1-m.loss[n])

    def R_balance_reservoir(self,model,a,t):
        """Reservoir temporal water balance."""
        m = model
        inf = m.U[a,t]
        ret = m.R[a,t]
        evp = m.V[a,t]*m.evap['incpt',a] + m.evap['slope']






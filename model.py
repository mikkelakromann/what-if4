from pyomo.environ import Objective, Constraint, Expression, Var, Set, Param
from pyomo.environ import NonNegativeReals
from pyomo.environ import SolverFactory, ConcreteModel, BuildAction
from pyomo.opt import SolverStatus, TerminationCondition


class whatif():
    """Singleton class for program instantiation"""

    def __init__(self):
        self.input = {}
        self.model = get_whatif_model(self.input.get_as_dict())
        self.model.solve()
        self.output = {}


def get_whatif_model(input_data):
    """Return fully initialised Pyomo model components: Sets, equations, expressions, variables, parameters."""

    m = ConcreteModel()

    sets = input_data['sets']
    R = m.R = Set(dimen=1, initialize=sets['R'])                    # Regions (eg. countries)
    N = m.N = Set(dimen=1, initialize=sets['N'])                    # River nodes (eg. up-, mid-, & downstream)
    T = m.T = Set(dimen=1, initialize=sets['T'])                    # Time tuples (eg. (y2020,m01,day) )
    C = m.C = Set(dimen=1, initialize=sets['C'])                    # Cultivation patterns (eg. wheat_cotton, ...)
    P = m.P = Set(dimen=1, initialise=sets['P'])                    # Phases of crop growth
    G = m.G = Set(dimen=1, initialize=sets['G'])                    # Goods (eg. elect, maize, wheat, sugar, ...)
    A = m.A = Set(dimen=1, initialize=sets['A'])                    # Assets (eg. turb1_100MW_rogun, irrig100ha, )

    subsets = input_data['subsets']
    # Time subsets
    TY = m.TY = Set(dimen=1, within=T, initialize=subsets['TY'])    # Time-years (e.g. y2020,mAll,hAll)
    TM = m.TM = Set(dimen=1, within=T, initialize=subsets['TM'])    # Time-months (e.g. y2020,m012,hAll)
    TH = m.TH = Set(dimen=1, within=T, initialize=subsets['TH'])    # Time-hours (e.g. y2020,m012,hDay)
    # Asset subsets
    AT = m.AT = Set(dimen=1, within=A, initialize=subsets['AT'])    # Assets for trade (transmission, transport)
    AE = m.AE = Set(dimen=1, within=A, initialize=subsets['AE'])    # Assets for electricity (thermal, solar, wind)
    AA = m.AA = Set(dimen=1, within=A, initialize=subsets['AA'])    # Assets for agriculture (rainfed+irig. fields)
    AH = m.AH = Set(dimen=1, within=A, initialize=subsets['AH'])    # Assets hydropower turbines
    AR = m.AR = Set(dimen=1, within=A, initialize=subsets['AR'])    # Assets reservoirs
    # Goods subsets
    GA = m.GA = Set(dimen=1, within=G, initialize=subsets['GA'])    # Goods from agriculture
    GE = m.GE = Set(dimen=1, within=G, initialize=subsets['GE'])    # Goods from energy = electricity
    # Nodes subset
    NE = m.NE = Set(dimen=1, within=G, initialize=subsets['NE'])    # Nodes with environmental constraints

    # Variables - capital formation
    m.K = Var(A,TY, within=NonNegativeReals)                        # Capital formation, investment in time-year
    # Variables - market trade of goods
    # TODO: Fix market variables to zero by good and time, e.g. only annual crops, hourly electricity
    m.S = Var(A,G,T, within=NonNegativeReals)                       # Supply of goods from asset at time
    m.D = Var(R,G,T, within=NonNegativeReals)                       # Demand in region for goods at time
    m.X = Var(A,G,T, within=NonNegativeReals)                       # Trade (export/import) of good by asset & time
    # Variables - water flows and use
    m.F = Var(N,TM, within=NonNegativeReals)                        # Outflow from node at time
    m.V = Var(A,TM, within=NonNegativeReals)                        # Volume of reservoir by reservoir asset & time
    m.U = Var(A,TM, within=NonNegativeReals)                        # Use of water (gross) by asset and time
    m.R = Var(A,TM, within=NonNegativeReals)                        # Return of water (gross) by asset and time
    # Variables - agriculture
    m.L = Var(A,C,TY, within=NonNegativeReals)                      # Land use by asset, crop mix and time
    m.W = Var(A,G,P,TY, within=NonNegativeReals)                    # Water used by plants and growth phase

    # Parameters
    m.grvity = 9.82                                                 # Gravity acceleration m/s^2
    para = input_data['para']
    m.ast = Param(m.capa_asst,A, initialize=para['asst'])           # Asset properties
    m.discnt = Param(TY, initialize=para['discnt'])
    m.trdlos = Param(AT,G, initialize=para['trdlos'])               # Loss trading in goods between markets
    m.demicp = Param(R,G, initialize=para['demicp'])
    m.demslp = Param(R,G, initialize=para['demslp'])
    m.demmax = Param(R,G, initialize=para['demmax'])
    m.wtrsup = Param(N,TM, initialize=para['wtrsup'])               # Water supply by node and month
    m.wtrlos = Param(N,TM, initialize=para['wtrlos'])               # Water loss through node by month
    m.flwmax = Param(N,TM, initialize=para['flwmax'])               # Environmental flows, maximum limit
    m.flwmin = Param(N,TM, initialize=para['flwmin'])               # Environmental flows, mniimum requirement
    m.elemax = Param(AE,TH, initialize=para['elemax'])              # Electricity availability, shr of max capacity
    m.crp_kY = Param(GA,P, initialise=para['crp_kY'])               # Crop yield water response by phase
    m.crpeto = Param(AA,GA,P, initialise=para['crpeto'])            # Optimal evapotranspiration mm/phase
    m.crpyld = Param(AA,GA, initialize=para['crpyld'])              # Optimal crop yield
    m.precip = Param(AA,TM, initialize=para['precip'])              # Precipitation on crops
    m.culcst = Param(AA,C, initialize=para['culcst'])
    m.phwgth = Param(GA,P,TM, initialize=para['phwgth'])            #

    rules = whatif_rules()
    # Expressions - usuable in model constraints and post-processing results
    m.CVAL = Expression(R,G,TY, rule=rules.expr_consumer_value)
    m.SUPL = Expression(R,G,TY, rule=rules.expr_market_supply)
    m.EXPT = Expression(R,G,TY, rule=rules.expr_market_export)
    m.IMPT = Expression(R,G,TY, rule=rules.expr_market_import)
    m.EVAP = Expression(AR,TM, rule=rules.expr_res_evapor)
    m.HEAD = Expression(AR,TM, rule=rules.expr_hpp_tbhead)
    # Objective
    m.objective = Objective(rule=rules.objective)
    # Constraints
    m.C_ast_capmax = Constraint(A,TY, rule=rules.cstr_ast_capmax)
    m.C_crp_balnce = Constraint(R,GA,TY, rule=rules.cstr_mkt_balnce)
    m.C_ele_balnce = Constraint(R,GE,TH, rule=rules.cstr_mkt_balnce)
    m.C_wtr_balnce = Constraint(N,TM, rule=rules.cstr_wtr_balnce)
    m.C_wtr_flwmax = Constraint(NE,TM, rule=rules.cstr_wtr_flwmax)
    m.C_wtr_flwmin = Constraint(NE,TM, rule=rules.cstr_wtr_flwmin)
    m.C_res_temprl = Constraint(AR,TM, rule=rules.cstr_res_temprl)
    m.C_res_evapor = Constraint(AR,TM, rule=rules.cstr_res_evapor)
    m.C_hpp_supply = Constraint(AH,GE,TM, rule=rules.cstr_hpp_supply)
    m.C_hpp_elemax = Constraint(AH,GE,TH, rule=rules.cstr_hpp_elemax)
    m.C_ele_supply = Constraint(AE,GE,TH, rule=rules.cstr_ele_supply)
    m.C_crp_supply = Constraint(AA,GA,TY, rule=rules.cstr_crp_supply)
    m.C_crp_wtruse = Constraint(AA,TM, rule=rules.cstr_crp_wtruse)
    m.C_crp_lnduse = Constraint(AA,TY, rule=rules.cstr_crp_lnduse)
    # Return Pyomo model object
    return m


class whatif_rules():
    """Container class for the WHAT-IF rules for constraints and expressions."""

    def objective(self,model):
        """Objective function."""
        m = model
        value = sum(m.discnt[ty]*m.CVAL[r,g,ty] for r in m.R for g in m.G for ty in m.TY)
        capex = sum(m.discnt[ty]*m.K[a,ty] for a in m.A for ty in m.TY)
        fopex = sum(m.discnt[ty]*m.L[aa,c,ty]*m.culcst[aa,c] for aa in m.AA for c in m.C for ty in m.TY)
        vopex = 0
        return value - capex - fopex - vopex

    def expr_consumer_value(self,model,r,g,ty):
        """Annual consumer valuation of electricity."""
        m = model
        slope = m.demslp[r,g]
        itcpt = m.demicp[r,g]
        # Return yearly value converted from hourly or monthly value
        if m.tfrq[g] == 'hourly':
            cval = sum((m.D[r,g,th]*slope + itcpt) * m.D[r,g,th] for th in m.T_T['thInTy',ty])
        elif m.tfrq[g] == 'monthly':
            cval = sum((m.D[r,g,tm]*slope + itcpt) * m.D[r,g,tm] for tm in m.T_T['tmInTy',ty])
        elif m.tfrq[g] == 'yearly':
            cval = (m.D[r,g,ty]*slope + itcpt) * m.D[r,g,ty]
        else:
            cval = 0
        return cval

    def cstr_ast_capmax(self,model,a,ty):
        """Maximum investment in assets."""
        model.K[a,ty] < model.ast['invMax',a]

    def cstr_mkt_balnce(self,model,r,g,t):
        """Goods market balance rule."""
        # Supply plus import equals export + demand
        return model.SUPL[r,g,t] + model.IMPT[r,g,t] == model.EXPT[r,g,t] + model.D[r,g,t]

    def expr_market_supply(self,model,r,g,t):
        """Market supply."""
        return sum(model.S[a,g,t] for a in model.A_R['S',g,r])

    def expr_market_export(self,model,r,g,t):
        """Market export."""
        return sum(model.X[a,g,t] for a in model.A_R['X',g,r])

    def expr_market_import(self,model,r,g,t):
        """Market import."""
        return sum(model.X[a,g,t]*(1-model.trdlos[a,g]) for a in model.A_R['I',g,r])

    def cstr_wtr_balnce(self,model,n,tm):
        """Node water balance."""
        m = model
        sup = m.wtrsup[n,tm]                                        # Water supply from inside node (rainfall/runoff)
        ups = sum(m.F[nu,tm] for nu in m.N_N[n])                    # Flow from upstream node by time
        use = sum(m.U[a,tm] for a in m.A_N['U',n])                  # Use of water by asset from node by time
        ret = sum(m.R[a,tm] for a in m.A_N['R',n])                  # Return of water by asset from node by time
        eff = 1-m.wtrlos[n]                                         # Efficiency: Water not lost along the river node
        return m.F[n,tm] == (ups + sup + ret - use)*eff             # River node loss is measured as share of inflow

    def cstr_res_temprl(self,model,ar,tm):            # Reservoir (+lakes+wetlands) are assets, not nodes
        """Reservoir temporal water balance."""
        m = model
        inf = m.U[ar,tm]                                            # Inflow is reservoir use from upstream node
        use = sum(m.U[a,tm] for a in m.A_A['frmRes',ar])            # Other assets' use of water from reservoir
        dis = m.R[ar,tm]                                            # Return is sum of turbine discharge and spills
        evp = m.EVAP[ar,tm]                                         # Evaporation (from expression below)
        nxt = sum(m.V[ar,t_nxt] for t_nxt in m.T_T['nxtM',tm])      # Next month-time segment
        return nxt == inf - dis - evp - use

    def expr_res_evapor(self,model,ar,tm):
        """Reservoir evaporation."""
        m = model
        # For now reservoir evaporation is linear F(volume) (TODO: Better!)
        area = m.V[ar,tm]*m.ast['resIcp',ar] + m.ast['resSlp',ar]
        return area*m.ast['resEvp']

    def cstr_wtr_flwmax(self,model,n_emax,tm):
        """Environmental maximum flows."""
        return model.F(n_emax,tm) < model.flwmax(n_emax,tm)

    def cstr_wtr_flwmin(self,model,n_emax,tm):
        """Environmental maximum flows."""
        return model.F(n_emax,tm) > model.flwmin(n_emax,tm)

    def cstr_hpp_wtrmax(self,model,ah,ge,tm):
        """Bernouilli constraint for turbines."""
        return 0

    def cstr_hpp_supply(self,model,ah,ge,tm):
        """Supply of good electricity from hydro power."""
        m = model
        # Return water to next node act as discharge through turbines
        # Spills are modelled as outflows from upstream node
        poten_energy = m.R[ah,tm] * m.HEAD[ah,tm] * m.gravity
        elect_energy = sum(m.S[ah,ge,th] for th in m.T_T['thInTm',tm]) / m.ast['hppEff',ah]
        return poten_energy == elect_energy

    def expr_hpp_tbhead(self,model,ah,tm):
        """Turbine head."""
        m = model
        if m.ast['hppRoR',ah]:
            head = m.ast['hppHea',ah]
        elif m.ast['hppRes',ah]:
            volu = sum(m.V[ar,tm] for ar in m.A_A['resHpp',ah])
            head = volu*m.ast['hppSlp',ah] + m.ast['hppIcp',ah]
        else:
            head = 0
        return head

    def cstr_hpp_elemax(self,model,ah,ge,th):
        """Maximum electricity production from HPP turbine."""
        m = model
        return m.S[ah,ge,th] < m.ast['hpMaxE',ah]*sum(m.K[ah,ty] for ty in m.T_T['tyInTh',th])

    def cstr_ele_supply(self,model,ae,ge,th):
        """Other electricity supply external availability limitations."""
        capacity = sum(model.K[ae,ty] for ty in model.T_T['tyInTh',th])
        model.S[ae,ge,th] < model.elemax[ae,th]*capacity

    def cstr_crp_supply(self,model,aa,ga,ty):
        """"Supply of agricultural crops including Yield-Water-Response."""
        m = model
        lndUse = sum(m.L[aa,c,ty] for c in m.C_G[ga])
        # Land deficit: The loss of crops in terms of ha land equivalent lost to under-evapotranspiration
        # W (unit m3) divided with crop_eto (unit m water) is equal to land use with optimal evapotranspiration
        # in which case land deficit is zero. If W is less than optimal, land deficit increases
        # Note: crop_Ky should be weighted with phase length as a share of total crop growth season
        #       i.e. if phase is 25% of time and unweighted kY is 0.80, the weighted kY is 0.20
        lndDef = sum(m.crp_kY[ga,p]*(lndUse - m.W[aa,ga,p,ty]/m.crpeto[aa,ga,p]) for p in m.P)
        # Supply is optimal yield x land use minus land deficit
        return m.S[aa,ga,ty] == m.crpyld[aa,ga] * (lndUse - lndDef)

    def cstr_crp_wtruse(self,model,aa,tm):
        """Crop water use including precipitation."""
        m = model
        use = sum(m.W[aa,ga,p,ty]*m.crop_pwght[ga,p,tm] for ga in m.GA for p in m.P for ty in m.T_T['tyInTm',tm])
        # If precipitation exceeds crop water need, abstraction of river water can be zero
        return m.U(aa,tm) + m.precip[aa,tm] > use

    def cstr_crp_lnduse(self,model,aa,ty):
        """Land use is limited by available land assets."""
        return sum(model.L[aa,c,ty] for c in model.C) < model.K[aa,ty]

    # end of file

"""
Custom Combine Physics Model for Charged Higgs Search: t -> H+ b, H+ -> W+ A(-> bb).
Implements non-linear decay branching ratio relations:
  POI: BR = B(t -> H+ b)
  tt_scale: (1 - BR) * (1 - BR)
  signal_scale: 2 * sigma_tt * (2 * 3 * 0.1086 * 0.6741 + 3 * 0.1086 * 3 * 0.1086) * (1 - BR) * BR
"""

from HiggsAnalysis.CombinedLimit.PhysicsModel import PhysicsModel


class ChargedHiggsBR(PhysicsModel):
    """Physics model parameterizing SM ttbar and Charged Higgs in terms of BR(t -> H+ b)."""

    def __init__(self):
        super(ChargedHiggsBR, self).__init__()
        self.sigma_tt = 923.6
        # W decay branching factors:
        # W -> l nu (per flavor): 0.1086, W -> hadrons: 0.6741
        # semileptonic: 2 * (3 * 0.1086) * 0.6741
        # dileptonic: (3 * 0.1086) * (3 * 0.1086)
        self.lep_filter_factor = (
            2.0 * 3.0 * 0.1086 * 0.6741 + (3.0 * 0.1086) * (3.0 * 0.1086)
        )
        self.signal_prefactor = 2.0 * self.sigma_tt * self.lep_filter_factor

    def setPhysicsOptions(self, physOptions):
        pass

    def doParametersOfInterest(self):
        # Define POI: BR = B(t -> H+ b) with range [0.0, 1.0]
        self.modelBuilder.doVar("BR[0.001,0.0,1.0]")
        self.modelBuilder.doSet("POI", "BR")

        # TTbar background scaling: (1 - BR) * (1 - BR)
        self.modelBuilder.factory_('expr::tt_scaling("(1.0 - @0)*(1.0 - @0)", BR)')

        # Signal scaling: 2 * sigma_tt * lep_filter_factor * (1 - BR) * BR
        self.modelBuilder.factory_(
            f'expr::signal_scaling("{self.signal_prefactor:.8f} * (1.0 - @0) * @0", BR)'
        )

    def getYieldScale(self, bin, process):
        if process == "signal":
            return "signal_scaling"
        elif process.startswith("TT_"):
            return "tt_scaling"
        else:
            return 1


chargedHiggsBR = ChargedHiggsBR()

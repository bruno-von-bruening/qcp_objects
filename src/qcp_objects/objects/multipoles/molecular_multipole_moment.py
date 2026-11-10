from . import *
from qcp_objects.crude_objects.multipoles.multipole_tensor import molecular_multipoles_base


class allowed_centers( str, Enum):
    COM = 'COM' # Center of Mass
    CONC = 'CONC' # Center of Nuclear Charge

class MolecularMultipoleMoments(molecular_multipoles_base):
    """Molecular multipole moment object.

    Attributes:
        center (str): Center of the multipole moment.
        representation (allowed_representation): Representation of the multipole moment, either 'S' for spherical or 'C' for cartesian.
    """
    __allowed_centers__ = allowed_centers
    expansion_center: Tuple[float,float,float]
    type_of_center: allowed_centers|None=None

    # Maybe I need internal rzepresentation



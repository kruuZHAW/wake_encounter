import sys
from pathlib import Path  # if you haven't already done so

# Add the root of this module to the path. This is necessary for this module to
# be imported from outside the package, see:
# https://stackoverflow.com/questions/16981921/relative-imports-in-python-3
file = Path(__file__).resolve()
parent, root = file.parent, file.parents[1]
sys.path.append(str(root))

from . import aircraft, meteo, wake

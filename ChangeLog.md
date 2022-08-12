# Summer 2022 Change Log
---
## Goals:
- ~~Move BSR & DBSR from version 3 libraries to version 4 libraries (libs_04)~~
- ~~Move executables from version 3 to version 4 (if applicable)~~
- Test all executables and utilities to verify compabitability (and changes if applicable)
- Refactor MPI implementation
- Increase performance by further implementing CUDA
- Boost uniformity of code structure
- ~~Reorganize repository to combine serial/MPI versions~~
---
## Change Log:
- Restructured respository so code is at same directory depth and similar patterning
- Removed version numbers from library and executables
- Added MPI functionality for BSR executables: CONF, BREIT
- Started implementing MPI for BSR executables: MAT, RECOUP, PHOT, HD
- Tested all serial BSR and DBSR latest version executables with version 4 library (exceptions: BSR_RECOUP, DBSR_jj2ls)
- Tested 34/46 utilities
- **Changed call structure of BSR_MULT: `bsr_mult initial=<initial state> final=<final state> AA=<interaction>`**
- Added optional CUDA library calls instead of LAPACK library calls (if available)
- Added CUDA detection in CMAKE
---
See commit messages and/or Trello board for more details
---
---

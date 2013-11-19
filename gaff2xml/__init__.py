#!/usr/bin/env python

"""
gaff2xml, a utility for using GAFF files in OpenMM

COPYRIGHT AND LICENSE

@author John D. Chodera <jchodera@gmail.com>
@author Kyle A. Beauchamp <kyleabeauchamp@gmail.com>

All code in this repository is released under the GNU General Public License.

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 
You should have received a copy of the GNU General Public License along with
this program.  If not, see <http://www.gnu.org/licenses/>.

"""

from gaff2xml import amber_parser, gafftools, system_checker

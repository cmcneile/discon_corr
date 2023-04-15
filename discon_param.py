#
#  QED parameters
#

import gvar as gv

w0 = gv.gvar('0.1715(9)')
hbarc = 0.197326968

w0overa_v_coarse = gv.gvar('1.1367(5)')  # very coarse ensemble    
w0overa_coarse  = gv.gvar('1.4149(6)')  # coarse ensemble
w0overa_fine = gv.gvar('1.9518(7)')  # coarse ensemble ???

#  Taste singlet ZV
ZV_v_coarse = gv.gvar('0.95932(18)') # 2005.01845 
ZV_coarse = gv.gvar('0.97255(22)') #  table X in 2005.0184
ZV_fine = gv.gvar('0.98445(11)') #  table X in 2005.01845

ZVqed_v_coarse = gv.gvar('0.999544(14)')*ZV_v_coarse # 2005.01845  
ZVqed_coarse = gv.gvar('0.999631(24)')*ZV_coarse #  table X in 2005.01845
ZVqed_fine = gv.gvar('0.999756(32)')*ZV_fine #  table X in 2005.01845    

#
#  taste singlet operator
#

#
#HPQCD, D. Hatton, C. T. H. Davies, G. P. Lepage, and A. T. Lytle, 
#Phys. Rev. D 100, 114513 (2019), arXiv:1909.00756.

Zvc = gv.gvar('0.93516(16)')/gv.gvar('0.820192(14)')
Zc = gv.gvar('0.94966(20)')/gv.gvar('0.834613(14)')
Zf = gv.gvar('0.96695(11)')/gv.gvar('0.852477(9)')

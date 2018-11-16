""""
  Wrapper class to store the loop operators for a single
  configuration.
"""

import copy

import discutil
##from lattice import *
from my_global import *
import loop_routines_stoch
import loop_routines_eigen

import sys

class loop_store:


    def __init__(self,wot, nt,nsample,tag_hhh, tag_eigen):
       self.wot = wot
       self.ns = 0
       self.nt = nt
       self.tag = tag_hhh 
       if wot == Vx_op  or  wot == Vy_op or  wot == Vz_op  :
          self.loop_stoch_l_1 = loop_routines_stoch.single_loop(nt, nsample,tag_hhh)
          self.loop_eigen_l_1 = loop_routines_eigen.single_loop_eigen(nt, tag_eigen)

          self.loop_stoch_s_1 = loop_routines_stoch.single_loop(nt, nsample,tag_hhh)
          self.loop_eigen_s_1 = loop_routines_eigen.single_loop_eigen(nt, tag_eigen)
       elif wot == Vxyz_op :
          self.loop_stoch_l_1 = loop_routines_stoch.single_loop(nt, nsample,tag_hhh)
          self.loop_eigen_l_1 = loop_routines_eigen.single_loop_eigen(nt, tag_eigen)
          self.loop_stoch_s_1 = loop_routines_stoch.single_loop(nt, nsample,tag_hhh)
          self.loop_eigen_s_1 = loop_routines_eigen.single_loop_eigen(nt, tag_eigen)

          self.loop_stoch_l_2 = loop_routines_stoch.single_loop(nt, nsample,"DIS_LOOP_Y_hpqcd")
          self.loop_eigen_l_2 = loop_routines_eigen.single_loop_eigen(nt, "EIGEN_DIS_LOOP_Y")
          self.loop_stoch_s_2 = loop_routines_stoch.single_loop(nt, nsample,"DIS_LOOP_Y_hpqcd")
          self.loop_eigen_s_2 = loop_routines_eigen.single_loop_eigen(nt, "EIGEN_DIS_LOOP_Y")

          self.loop_stoch_l_3 = loop_routines_stoch.single_loop(nt, nsample,"DIS_LOOP_Z_hpqcd")
          self.loop_eigen_l_3 = loop_routines_eigen.single_loop_eigen(nt,"EIGEN_DIS_LOOP_Z")
          self.loop_stoch_s_3 = loop_routines_stoch.single_loop(nt, nsample,"DIS_LOOP_Z_hpqcd")
          self.loop_eigen_s_3 = loop_routines_eigen.single_loop_eigen(nt, "EIGEN_DIS_LOOP_Z")
       else:
          print ("Error wot = " , wot)
          sys.exit(0)

    def loadin(self, filename, mass_light, mass_strange) :
       wot = self.wot

       if wot == Vx_op or  wot == Vy_op  or  wot == Vz_op   :
          self.loop_stoch_l_1.load_data_mass(filename,mass_light) 
          self.loop_eigen_l_1.load_data_vec_im_mass(filename,mass_light) 

          self.loop_stoch_s_1.load_data_mass(filename,mass_strange) 
          self.loop_eigen_s_1.load_data_vec_im_mass(filename,mass_strange) 
       elif wot == Vxyz_op :
          self.loop_stoch_l_1.load_data_mass(filename,mass_light) 
          self.loop_eigen_l_1.load_data_vec_im_mass(filename,mass_light) 

          self.loop_stoch_s_1.load_data_mass(filename,mass_strange) 
          self.loop_eigen_s_1.load_data_vec_im_mass(filename,mass_strange) 

          self.loop_stoch_l_2.load_data_mass(filename,mass_light) 
          self.loop_eigen_l_2.load_data_vec_im_mass(filename,mass_light) 

          self.loop_stoch_s_2.load_data_mass(filename,mass_strange) 
          self.loop_eigen_s_2.load_data_vec_im_mass(filename,mass_strange) 

          self.loop_stoch_l_3.load_data_mass(filename,mass_light) 
          self.loop_eigen_l_3.load_data_vec_im_mass(filename,mass_light) 

          self.loop_stoch_s_3.load_data_mass(filename,mass_strange) 
          self.loop_eigen_s_3.load_data_vec_im_mass(filename,mass_strange) 
       else:
          print ("Error wot = " , wot )
          sys.exit(0)




    def compute(self, STATns):
       wot = self.wot

       if wot == Vx_op  or  wot == Vy_op  or  wot == Vz_op   :
          self.loop_stoch_l_1.adjust_samples_compute(STATns)
          self.loop_stoch_l_1.compute()
          self.loop_eigen_l_1.compute()

          self.loop_stoch_s_1.adjust_samples_compute(STATns)
          self.loop_stoch_s_1.compute()
          self.loop_eigen_s_1.compute()
       elif wot == Vxyz_op :
          self.loop_stoch_l_1.adjust_samples_compute(STATns)
          self.loop_stoch_l_1.compute()
          self.loop_eigen_l_1.compute()

          self.loop_stoch_s_1.adjust_samples_compute(STATns)
          self.loop_stoch_s_1.compute()
          self.loop_eigen_s_1.compute()
          
          self.loop_stoch_l_2.adjust_samples_compute(STATns)
          self.loop_stoch_l_2.compute()
          self.loop_eigen_l_2.compute()

          self.loop_stoch_s_2.adjust_samples_compute(STATns)
          self.loop_stoch_s_2.compute()
          self.loop_eigen_s_2.compute()
          
          self.loop_stoch_l_3.adjust_samples_compute(STATns)
          self.loop_stoch_l_3.compute()
          self.loop_eigen_l_3.compute()

          self.loop_stoch_s_3.adjust_samples_compute(STATns)
          self.loop_stoch_s_3.compute()
          self.loop_eigen_s_3.compute()
          
       else:
          print ("Error wot = " , wot )
          sys.exit(0)




    def include_eigen_corr(self) :
       wot = self.wot

       if wot == Vx_op or  wot == Vy_op or  wot == Vz_op  :
          self.loop_stoch_l_1.include_eigen_corr(self.loop_eigen_l_1) 
          self.loop_stoch_s_1.include_eigen_corr(self.loop_eigen_s_1) 
       elif wot == Vxyz_op :
          self.loop_stoch_l_1.include_eigen_corr(self.loop_eigen_l_1) 
          self.loop_stoch_s_1.include_eigen_corr(self.loop_eigen_s_1) 

          self.loop_stoch_l_2.include_eigen_corr(self.loop_eigen_l_2) 
          self.loop_stoch_s_2.include_eigen_corr(self.loop_eigen_s_2) 

          self.loop_stoch_l_3.include_eigen_corr(self.loop_eigen_l_3) 
          self.loop_stoch_s_3.include_eigen_corr(self.loop_eigen_s_3) 
       else:
          print ("Error wot = " , wot )
          sys.exit(0)




    def include_eigen_corr_fact(self, fact) :
       wot = self.wot

       if wot == Vx_op or  wot == Vy_op or  wot == Vz_op  :
          self.loop_stoch_l_1.include_eigen_corr_fact(self.loop_eigen_l_1,fact) 
          self.loop_stoch_s_1.include_eigen_corr_fact(self.loop_eigen_s_1,fact) 
       elif wot == Vxyz_op :
          self.loop_stoch_l_1.include_eigen_corr_fact(self.loop_eigen_l_1,fact) 
          self.loop_stoch_s_1.include_eigen_corr_fact(self.loop_eigen_s_1,fact) 

          self.loop_stoch_l_2.include_eigen_corr_fact(self.loop_eigen_l_2,fact) 
          self.loop_stoch_s_2.include_eigen_corr_fact(self.loop_eigen_s_2,fact) 

          self.loop_stoch_l_3.include_eigen_corr_fact(self.loop_eigen_l_3,fact) 
          self.loop_stoch_s_3.include_eigen_corr_fact(self.loop_eigen_s_3,fact) 
       else:
          print ("Error wot = " , wot )
          sys.exit(0)



    def replace_with_eigen_corr(self) :
       wot = self.wot

       if wot == Vx_op or  wot == Vy_op or  wot == Vz_op  :
          self.loop_stoch_l_1.replace_with_eigen_corr(self.loop_eigen_l_1) 
          self.loop_stoch_s_1.replace_with_eigen_corr(self.loop_eigen_s_1) 
       else:
          print ("Error wot = " , wot )
          sys.exit(0)



    def create_light_strange(self):

       if self.wot == Vx_op or  self.wot == Vy_op or  self.wot == Vz_op  :
           self.loop_stoch_1 = copy.deepcopy(self.loop_stoch_l_1)
           self.loop_stoch_1.subtract(self.loop_stoch_s_1)
       elif self.wot == Vxyz_op : 
           self.loop_stoch_1 = copy.deepcopy(self.loop_stoch_l_1)
           self.loop_stoch_1.subtract(self.loop_stoch_s_1)

           self.loop_stoch_2 = copy.deepcopy(self.loop_stoch_l_2)
           self.loop_stoch_2.subtract(self.loop_stoch_s_2)

           self.loop_stoch_3 = copy.deepcopy(self.loop_stoch_l_3)
           self.loop_stoch_3.subtract(self.loop_stoch_s_3)

       else:
          print ("Error wot = " , wot )
          sys.exit(0)



    def dump_loop(self):

       if self.wot == Vx_op or  self.wot == Vy_op or  self.wot == Vz_op  :
           print ("Light stochastic loops " )
           self.loop_stoch_l_1.dump_time_avloop()

           print ( "Strange stochastic loops " )
           self.loop_stoch_s_1.dump_time_avloop()
       else:
          print ("dump_loop  Error wot = " , wot )
          sys.exit(0)



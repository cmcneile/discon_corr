#  Utility functions
#
#

from my_global import *

def set_tags_twisted(wot) :

   if wot == Vx_op :
      tag_hhh ="DIS_LOOP_X_TWISTED"
      tag_eigen="EIGEN_DIS_LOOP_X"
      title = "Vector_X"
      out_tag = "summary_v1" 
   elif wot == unit_op :
      tag_hhh = "DIS_LOOP_one_TWISTED"
      tag_eigen="EIGEN_DIS_LOOP_unit"
      title = "Gamma_one"
      out_tag = "summary_unit" 
   elif wot == Vxyz_op :
      tag_hhh ="DIS_LOOP_X_TWISTED"
      tag_eigen="EIGEN_DIS_LOOP_X"
      title = "Vector_X + Vector_Y + Vector_Z"
      out_tag = "summary_v1v2v3" 
   elif wot == Vy_op  :
      tag_hhh ="DIS_LOOP_Y_TWISTED"
      tag_eigen="EIGEN_DIS_LOOP_Y"
      title = "Vector_Y"
      out_tag = "summary_v2" 
   elif wot == Vz_op :
      tag_hhh ="DIS_LOOP_Z_TWISTED"
      tag_eigen="EIGEN_DIS_LOOP_Z"
      title = "Vector_Z"
      out_tag = "summary_v3" 
   else :
      print ("ERROR Analysis wot = " , wot , " out of range " )
      sys.exit(0)

   return tag_hhh, tag_eigen, title , out_tag 


def set_tags_twisted_volume(wot) :

   if wot == Vx_op :
      tag_hhh ="DIS_LOOP_X_Vxyzt_TWISTED"
      tag_eigen="EIGEN_DIS_LOOP_X"
      title = "Vector_X_volume"
      out_tag = "summary_v1" 
   elif wot == unit_op :
      tag_hhh = "DIS_LOOP_one_Vxyzt_TWISTED"
      tag_eigen="EIGEN_DIS_LOOP_unit"
      title = "Gamma_one_volume"
      out_tag = "summary_unit" 
   elif wot == Vxyz_op :
      tag_hhh ="DIS_LOOP_X_Vxyzt_TWISTED"
      tag_eigen="EIGEN_DIS_LOOP_X"
      title = "Vector_X + Vector_Y + Vector_Z"
      out_tag = "summary_v1v2v3" 
   elif wot == Vy_op  :
      tag_hhh ="DIS_LOOP_Y_Vxyzt_TWISTED"
      tag_eigen="EIGEN_DIS_LOOP_Y"
      title = "Vector_Y_volume"
      out_tag = "summary_v2" 
   elif wot == Vz_op :
      tag_hhh ="DIS_LOOP_Z_Vxyzt_TWISTED"
      tag_eigen="EIGEN_DIS_LOOP_Z"
      title = "Vector_Z_volume"
      out_tag = "summary_v3" 
   else :
      print ("ERROR Analysis wot = " , wot , " out of range " )
      sys.exit(0)

   return tag_hhh, tag_eigen, title , out_tag 






def set_tags(wot) :

   if wot == Vx_op  :
      tag_hhh ="DIS_LOOP_X_hpqcd"
      tag_eigen="EIGEN_DIS_LOOP_X"
      title = "Vector_X"
      out_tag = "summary_v1"
   elif wot == 1 :
      tag_hhh = "DIS_LOOP_one_hpqcd"
      tag_eigen="EIGEN_DIS_LOOP_unit"
      title = "Gamma_one"
      out_tag = "summary_unit"
   elif wot == Vxyz_op  :
      tag_hhh ="DIS_LOOP_X_hpqcd"
      tag_eigen="EIGEN_DIS_LOOP_X"
      title = "Vector_X + Vector_Y + Vector_Z"
      out_tag = "summary_v1v2v3"
   elif wot == Vy_op  :
      tag_hhh ="DIS_LOOP_Y_hpqcd"
      tag_eigen="EIGEN_DIS_LOOP_Y"
      title = "Vector_Y"
      out_tag = "summary_v2"
   elif wot == Vz_op :
      tag_hhh ="DIS_LOOP_Z_hpqcd"
      tag_eigen="EIGEN_DIS_LOOP_Z"
      title = "Vector_Z"
      out_tag = "summary_v3"
   else :
      print ("ERROR Analysis wot = " , wot , " out of range "  )
      sys.exit(0)
       

   return tag_hhh, tag_eigen, title , out_tag 

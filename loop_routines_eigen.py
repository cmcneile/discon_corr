
#
#

import re
import math
import numpy as np
import sys


#
#



#
#  Class for reading the  eigenvalue contribution to disconnected diagrams
#
#


class single_loop_eigen:
   def __init__(self,nt,tag):
       self.nt = nt 
       self.tag = tag
       self.corr = np.zeros( (nt)  ) 


#  load data for unit loop
   def load_data_unit(self,filename):
       self.load_data(filename,5,4)

#  load data for real parr of vector
   def load_data_vec_re(self,filename):
       self.load_data(filename,5,4)

#  load data for imaginary parr of vector
   def load_data_vec_im(self,filename):
       self.load_data(filename,6,4)

   def load_data_vec_im_mass(self,filename,mass_in):
       self.load_data_mass(filename,6,4,mass_in)

#  load data for imaginary parr of vector
   def load_data_vec_im_inc(self,filename,tag):
       self.load_data_inc(filename,6,4,tag)

   def load_data_vec_re_inc(self,filename,tag):
       self.load_data_inc(filename,5,4,tag)

   def load_data_vec_im_mass_inc(self,filename,tag, mass_in):
       self.load_data_mass_inc(filename,6,4,tag, mass_in)


#    Load the data
#
#    ic correlator row
#    it time row

   def load_data(self,filename,ic,it):
      tag = self.tag
      print ("Reading eigen-correlators from " , filename , " searching for " , tag )

      f = open(filename, 'rU')
      text = f.readlines()


      count = 0 
      
      for line in text:
         lll = line.rstrip()
         if re.search(tag,  lll):
            tmp = lll.split()
#            print "DEBUG " , tmp[it], tmp[ic]
            ccc = float(tmp[ic])

            tt = int(tmp[it])

            self.corr[tt] = ccc

      f.close()



#    Load the data
#
#    ic correlator row
#    it time row

   def load_data_mass(self,filename,ic,it,mass_in):
      tag = self.tag
      print ("Reading eigen-correlators from " , filename , " searching for " , tag, " mass= " , mass_in)

      f = open(filename, 'rU')
      text = f.readlines()


      count = 0 
      
      for line in text:
         lll = line.rstrip()
         if re.search(tag,  lll):
            tmp = lll.split()
            mass = float(tmp[2] )
            if mass == mass_in :  
               ccc = float(tmp[ic])
               tt = int(tmp[it])
               self.corr[tt] = ccc

      f.close()


   def load_data_inc(self,filename,ic,it,tag):
      print ("Reading (inc) eigen-correlators from " , filename , " searching for " , tag)

      f = open(filename, 'rU')
      text = f.readlines()


      count = 0 
      
      for line in text:
         lll = line.rstrip()
         if re.search(tag,  lll):
#            print lll
            tmp = lll.split()
#            print "DEBUG " , tmp[it], tmp[ic]
            ccc = float(tmp[ic])

            tt = int(tmp[it])

#            print tt , ii , ccc
            self.corr[tt] += ccc

      f.close()




   def load_data_mass_inc(self,filename,ic,it,tag, mass_in):
      print ("Reading (inc) eigen-correlators from " , filename , " searching for " , tag)

      f = open(filename, 'rU')
      text = f.readlines()


      count = 0 
      
      for line in text:
         lll = line.rstrip()
         if re.search(tag,  lll):
            tmp = lll.split()
            mass = float(tmp[2] )
            if mass == mass_in :  
              ccc = float(tmp[ic])
              tt = int(tmp[it])
              self.corr[tt] += ccc

      f.close()




##

   def dump_time_avcorr(self):

      print ("Loop (eigenvalues) operators for " , self.tag )
      for tt in range(0,self.nt):
         print (tt , "%e" % self.corr[tt] )

   def apply_norm(self,norm):

      print (self.tag, "normalizing eigen correlators with norm = " , norm   )
      for tt in range(0,self.nt):
            self.corr[tt] = self.corr[tt] / norm


   def compute(self):
      tot = 0 
      for tt in range(0,self.nt):
          tot = tot + self.corr[tt]
      self.corr_t     = tot


   def dump_time_sum(self):
         cc_sum = self.corr_t
         print ( "INFO TSUM_EIGEN", self.tag , "%f" % cc_sum  )



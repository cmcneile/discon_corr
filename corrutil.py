"""
 Class to store averaged correlators and their errors
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import jackknife
import sys
import h5py
import math

class corr_sum:
   """
   Class to store the summary of correlators
   """

   def __init__(self,nt,ns,tend, tag):

       self.ns = ns
       self.nt = nt
       self.tend = tend
       self.tag = tag 

       self.dis_corr_err  = np.zeros( (self.tend)  )
       self.fdis_corr = np.zeros( (self.tend)  )
       self.tt = np.zeros( (self.tend)  )

       self.all_corr = np.zeros( (self.tend, self.ns)  )


   def norm(self,nrm) :
       """
         Normalize the correlators
       """
       for t in range(0,self.tend) :
          self.fdis_corr[t]     /= nrm
          self.dis_corr_err[t] /= nrm
          for ii in range(0,self.ns) :
              self.all_corr[t][ii]  /= nrm


   def write_dis_loop(self):
      """
      Write out the disconnected loops to standard output

      """

      print ("Disconnected loop for " , self.tag)
      for t in range(0,self.tend):
          print (t , "%e" % self.floop[t])


   def loadin(self,corr,corr_err):
      """
       Load the correlators into the class
      """

      for t in range(0,self.tend):
          self.tt[t] = t
          self.dis_corr_err[t] = corr_err[t]
          self.fdis_corr[t] =   corr[t] 



   def write_discorr(self, verbose = True):
      """
      Write out the disconnected correlators with errors to standard output
      """

      if verbose :
        print ("Disconnected correlators for " , self.tag)
        print ("nt = " , self.nt )
        print ("nt = " , self.tend)
        print ("ns = " , self.ns)

      for t in range(0,self.tend):
          if verbose :
             print ("SINFO", t , "%e" % self.fdis_corr[t]   , "%e" % self.dis_corr_err[t] )
          else:
             print ( t , "%e" % self.fdis_corr[t]   , "%e" % self.dis_corr_err[t] )



   def check_discorr(self):
      """
      Check that the correlators averaged over configurations equals the full sample corre;lators
      """

      print ("Disconnected correlators for " , self.tag)
      for t in range(0,self.tend):
          corr_f = 0.0 
          for iss in range(0, self.ns) :
            corr_f = corr_f + self.all_corr[t][iss]
          corr_f /= self.ns

          ans = self.fdis_corr[t] - corr_f 
          print ("CHECK", t , "%e" % ans   )


   def plot_data(self,l_fmt,  l_label, t_off=0):
       """
          Plot the data
       """

       if t_off == 0 :
         plt.errorbar(self.tt, self.fdis_corr, self.dis_corr_err ,  fmt= l_fmt , label = l_label)
       else: 
          tmp  = np.zeros( (self.tend)  )
          for t in range(0, self.tend ) :
#              print ("DEBUG t= " , t )
              tmp[t] = self.tt[t] +  t_off
          plt.errorbar(tmp, self.fdis_corr, self.dis_corr_err ,  fmt= l_fmt , label = l_label)


   def plot_data_norm(self,l_fmt,  l_label, nrm,  t_off=0 ):
       """
          Plot the data
       """

       tmp  = np.zeros( (self.tend)  )
       corr = np.zeros( (self.tend)  )
       corr_err = np.zeros( (self.tend)  )

       for t in range(0, self.tend ) :
          tmp[t] = self.tt[t] +  t_off
          corr[t] =  self.fdis_corr[t] / nrm
          corr_err[t] = self.dis_corr_err[t]  / nrm
#          print("DEBUG " , tmp[t] , corr[t] , corr_err[t] ) 


       plt.errorbar(tmp, corr , corr_err  ,  fmt= l_fmt , label = l_label)
#       plt.errorbar(self.tt , corr , corr_err  ,  fmt= l_fmt , label = l_label)

#       for tt in tmp:
#          print ("debug " , tt)

       return corr


   def plot_signalTOnoise(self,l_fmt,  l_label,   t_off=0 ):
       """
          Plot the data
       """

       tmp  = np.zeros( (self.tend)  )
       corr = np.zeros( (self.tend)  )

       for t in range(0, self.tend ) :
          tmp[t] = self.tt[t] +  t_off
          corr[t] =  self.dis_corr_err[t] / math.fabs(self.fdis_corr[t] )

##       plt.plot(tmp, corr ,   fmt= l_fmt , label = l_label)
       plt.plot(tmp, corr , l_fmt , label = l_label )
#       plt.errorbar(self.tt , corr , corr_err  ,  fmt= l_fmt , label = l_label)

#       for tt in tmp:
#          print ("debug " , tt)




   def plot_error(self,l_fmt,  l_label, nrm,  t_off=0 ):
       """
          Plot the data
       """

       tmp  = np.zeros( (self.tend)  )
       corr_err = np.zeros( (self.tend)  )

       for t in range(0, self.tend ) :
          tmp[t] = self.tt[t] +  t_off
          corr_err[t] = self.dis_corr_err[t]  / nrm
##       plt.plot(tmp, corr ,   fmt= l_fmt , label = l_label)
       plt.plot(tmp, corr_err , l_fmt , label = l_label )
#       plt.errorbar(self.tt , corr , corr_err  ,  fmt= l_fmt , label = l_label)

#       for tt in tmp:
#          print ("debug " , tt)



   def plot_a_mu(self,l_fmt,  l_label, fname, Vspace, Zv, tend= 20 , t_off=0 ):
       """
          Plot the correlator times weight
       """


       tmp  = np.zeros( (tend)  )
       corr = np.zeros( (tend)  )
       corr_err = np.zeros( (tend)  )

       print ("Plot a_mu^{HVP LO} for " , self.tag)

       print ("Reloading the coefficients from " , fname )
       if not os.path.isfile(fname) :
          print ("Can not find " , fname)
          sys.exit(1)

       print ("SINFO Norm = (taste =4 ) * (spatial volume) * (charge = 1/9)   "  )
       print ("SINFO spatial volume = " , Vspace)
       print ("SINFO T = " , self.tend)
       print ("SINFO ZV = " , Zv)

       norm = Vspace * 4 * 9  / (Zv*Zv)
       
       cc = np.load(fname) 
       t = 0 

       if True :
          for c in cc: 
             print ("w(" , t , ") = " , c )
             t += 1 

       if len(cc) < tend  :
         print ("Not enough coefficients" )
         sys.exit(1)


       for t in range(0, tend ) :
          tmp[t] = self.tt[t] +  t_off
          corr[t]     = cc[t] * self.fdis_corr[t] / norm
          corr_err[t] = cc[t] * self.dis_corr_err[t]   / norm


       plt.errorbar(tmp, corr , corr_err  ,  fmt= l_fmt , label = l_label)



   def estimate_amu(self, fname, Vspace, t_star, Zv, Verbose=False):
      """
      Estimate the a_\mu^{HVP (LO) DISC} using the UKQCD/RBC method  (1610.04603).
      """

      print ("Estimate a_mu^{HVP LO}  for " , self.tag, " using t*= " , t_star )

      print ("Reloading the coefficients from " , fname )
      if not os.path.isfile(fname) :
         print ("Can not find " , fname)
         sys.exit(1)

      print ("SINFO Norm = (taste =4 ) * (spatial volume) * (charge = 1/9)   "  )
      print ("SINFO spatial volume = " , Vspace)
      print ("SINFO T = " , self.tend)
      print ("SINFO ZV = " , Zv)

      norm = Vspace * 4 * 9  / (Zv*Zv)

      cc = np.load(fname) 
      t = 0 

      if Verbose :
        for c in cc: 
           print ("w(" , t , ") = " , c )
           t += 1 

      if len(cc) < t_star :
         print ("Not enough coefficients" )
         sys.exit(1)

      ans = 0.0 
      for t in range(0,t_star ):
         ans +=  self.fdis_corr[t] * cc[t]          

      ans /= norm   


      #  jackknife analysis
      j_amu = np.zeros( (self.ns)  )
      for i_ss in range(0,self.ns) :
         j_amu[i_ss] = 0.0 
         for t in range(0,t_star):
            for j_ss in range(0,self.ns) :
                if i_ss != j_ss :
                   j_amu[i_ss] += self.all_corr[t][j_ss]  * cc[t]

         j_amu[i_ss]  /= (self.ns - 1)

      ans_err = jackknife.jackknife(j_amu,self.ns )
      ans_err /= norm 

      print ("SINFO a_mu^{DIS HVP LO}" ,  "%e" % ans , " +/- %e " % ans_err  )

      print ("ANS_A_MU " , t_star ,  "%e" % ans , " +/- %e " % ans_err  )

      return ans , ans_err




   def plot_estimate_amu(self, fname, Vspace, t_star_end, Zv, l_fmt, l_label, Verbose=False, sss=0.0):
      """
      Estimate the a_\mu^{HVP (LO) DISC} using the UKQCD/RBC method  (1610.04603).
      """

      print ("PLOT Estimate a_mu^{HVP LO}  for " , self.tag )

      print ("Reloading the coefficients from " , fname )
      if not os.path.isfile(fname) :
         print ("Can not find " , fname)
         sys.exit(1)

      print ("SINFO Norm = (taste =4 ) * (spatial volume) * (charge = 1/9)   "  )
      print ("SINFO spatial volume = " , Vspace)
      print ("SINFO T = " , self.tend)
      print ("SINFO ZV = " , Zv)

      norm = Vspace * 4 * 9  / (Zv*Zv)

      cc = np.load(fname) 
      t = 0 

      if Verbose :
        for c in cc: 
           print ("w(" , t , ") = " , c )
           t += 1 

      if len(cc) < t_star_end :
         print ("Not enough coefficients" )
         sys.exit(1)

      tmp  = np.zeros( (t_star_end)  )
      corr = np.zeros( (t_star_end)  )
      corr_err = np.zeros( (t_star_end)  )


      for t_star in range(1, t_star_end ) :
         ans = 0.0 
         for t in range(0,t_star ):
            ans +=  self.fdis_corr[t] * cc[t]          
         ans /= norm   

         #  jackknife analysis
         j_amu = np.zeros( (self.ns)  )
         for i_ss in range(0,self.ns) :
            j_amu[i_ss] = 0.0 
            for t in range(0,t_star):
               for j_ss in range(0,self.ns) :
                if i_ss != j_ss :
                   j_amu[i_ss] += self.all_corr[t][j_ss]  * cc[t]

            j_amu[i_ss]  /= (self.ns - 1)

         ans_err = jackknife.jackknife(j_amu,self.ns )
         ans_err /= norm 

         print ("ANS_A_MU " , t_star ,  "%e" % ans , " +/- %e " % ans_err  )
         tmp[t_star] = t_star + sss
         corr[t_star] = ans
         corr_err[t_star] = ans_err
     
      plt.errorbar(tmp, corr , corr_err  ,  fmt= l_fmt , label = l_label)

#     https://scipy.github.io/old-wiki/pages/Cookbook/InputOutput.html
      f_amu = "amu_estimate.npz"
      print ("Writing a_mu(T) to " , f_amu)
      np.savez(f_amu, t=tmp,corr=corr, corr_err=corr_err) 



   def plot_diff_estimate_amu(self, sub, fname, Vspace, t_star_end, Zv, l_fmt, l_label, Verbose=False):
      """
      Estimate the a_\mu^{HVP (LO) DISC} using the UKQCD/RBC method  (1610.04603).
      """

      print ("PLOT Estimate difference a_mu^{HVP LO}  for " , self.tag )

      print ("Reloading the coefficients from " , fname )
      if not os.path.isfile(fname) :
         print ("Can not find " , fname)
         sys.exit(1)

      print ("SINFO Norm = (taste =4 ) * (spatial volume) * (charge = 1/9)   "  )
      print ("SINFO spatial volume = " , Vspace)
      print ("SINFO T = " , self.tend)
      print ("SINFO ZV = " , Zv)

      norm = Vspace * 4 * 9  / (Zv*Zv)

      cc = np.load(fname) 
      t = 0 

      if Verbose :
        for c in cc: 
           print ("w(" , t , ") = " , c )
           t += 1 

      if len(cc) < t_star_end :
         print ("Not enough coefficients" )
         sys.exit(1)

      tmp  = np.zeros( (t_star_end)  )
      corr = np.zeros( (t_star_end)  )
      corr_err = np.zeros( (t_star_end)  )


      for t_star in range(1, t_star_end ) :
         ans = 0.0 
         for t in range(0,t_star ):
            ans +=  (self.fdis_corr[t] - sub.fdis_corr[t] ) * cc[t]          
         ans /= norm   

         #  jackknife analysis
         j_amu = np.zeros( (self.ns)  )
         for i_ss in range(0,self.ns) :
            j_amu[i_ss] = 0.0 
            for t in range(0,t_star):
               for j_ss in range(0,self.ns) :
                if i_ss != j_ss :
                   j_amu[i_ss] += (self.all_corr[t][j_ss] - sub.all_corr[t][j_ss] ) * cc[t]

            j_amu[i_ss]  /= (self.ns - 1)

         ans_err = jackknife.jackknife(j_amu,self.ns )
         ans_err /= norm 

         print ("ANS_A_MU " , t_star ,  "%e" % ans , " +/- %e " % ans_err  )
         tmp[t_star] = t_star
         corr[t_star] = ans
         corr_err[t_star] = ans_err
     
      plt.errorbar(tmp, corr , corr_err  ,  fmt= l_fmt , label = l_label)

      f_amu = "diff_amu_estimate.npz"
      print ("Writing a_mu(T) to " , f_amu)
      np.savez(f_amu, t=tmp,corr=corr, corr_err=corr_err) 


#  string hack from
#  https://stackoverflow.com/questions/23220513/storing-a-list-of-strings-to-a-hdf5-dataset-from-python

   def  export_to_h5(self, fname, cfglist) :

       dt = h5py.special_dtype(vlen=str)
       with h5py.File(fname,"w") as hf:
           hf.create_dataset('discon_corr', data=self.all_corr) 
#           hf.create_dataset('lattice_number', data=cfglist, dtype="S10") 
           asciiList = [n.encode("ascii", "ignore") for n in cfglist]

           hf.create_dataset('lattice_number', (len(asciiList),1),'S10', asciiList)

       print ("Correlators exported to the file " , fname)


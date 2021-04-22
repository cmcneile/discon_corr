"""
  Routines to create disconnected correlatots from loop
  measurements.

  The class stores the loops for all configurations
  ansd then constructs the correlators.

"""
import numpy as np
import pickle
import os.path
import sys

import jackknife
import corrutil

class dis_corr:
   """
   Class to create disconnected correlators from loop operators.

   """

   def __init__(self,nt,ns,tend, tag):
       self.ns = 0
       self.nt = nt
       self.tend = tend
       self.tag = tag 
       self.loop = [] 


   def load(self, xx):
        self.loop.append(xx)
        self.ns = self.ns + 1

   def  compute_jack_loop(self ):


       print ("Computing jackknife loops " , self.tag)
       self.jloop = np.zeros( (self.nt, self.ns)  )
       for t in range(0,self.nt):
         for iss in range(0,self.ns) :
            self.jloop[t][iss] = 0.0 
            for js in range(0,self.ns) :
               if js != iss :
                 self.jloop[t][iss] = self.jloop[t][iss] +  self.loop[js].av_loop[t]
            self.jloop[t][iss] = self.jloop[t][iss] / (self.ns - 1) 

   def  dump_jack_loop(self ):

       print ("Dump jack loops")
       for t in range(0,self.nt):
         for iss in range(0,self.ns) :
            print (t, iss, self.jloop[t][iss] )



   def  compute_full_loop(self ):

       print ("Computing full sample loops " , self.tag)
       self.floop = np.zeros( (self.nt)  )
       for t in range(0,self.nt):
         self.floop[t] =  0.0
         for iss in range(0,self.ns) :
            self.floop[t] = self.floop[t]+  self.loop[iss].av_loop[t]
         self.floop[t] = self.floop[t] / (self.ns) 

   def dump_full_loop(self ):

       print ("Full loop ")
       for t in range(0,self.nt):
           print ("DISCORR_F_DUMP" , t , self.floop[t])
   
   def compute(self):
       """
          Compute the full and jacksamples of the loops.
          Use the loops to compute the disconnected correlators.

       """

       self.compute_full_loop()
       self.compute_jack_loop()

       self.compute_full_discorr() 
       self.compute_discorr_err()


   def compute_unbias(self):
       """
          Compute the full and jacksamples of the loops.
          Use the loops to compute the disconnected correlators.

       """

       self.compute_full_loop()
       self.compute_jack_loop()

       self.compute_full_discorr_unbias() 
       self.compute_discorr_err()

   def compute_vacsub(self):

       self.compute_full_loop()
       self.compute_jack_loop()

       self.compute_full_discorr_vacsub() 
       self.compute_discorr_vacsub_err()



   def compute_full_discorr_A(self):
       """
          Compute the full sample correlators by first averaging over loops 
       """

       self.fdis_corr = np.zeros( (self.tend)  )

       for t in range(0,self.tend):
         self.fdis_corr[t] = 0.0
         for tt in range(0,self.nt) :
           tinc = t + tt
           if( tinc >= self.nt):
             tinc = tinc - self.nt
           self.fdis_corr[t] = self.fdis_corr[t] + self.floop[tt]*self.floop[tinc]
         self.fdis_corr[t] /=  self.nt


   def compute_full_discorr(self):
       """
          Compute the correlators on each comfiguration and then average
       """

       self.fdis_corr = np.zeros( (self.tend)  )

       for iss in range(0, self.ns) :
         self.loop[iss].compute_corr(self.tend) 
##         print ("DEBUG corr[0] " , self.loop[iss].av_corr[0] , self.loop[iss].av_corr[1] )
         for t in range(0,self.tend):
              self.fdis_corr[t] = self.fdis_corr[t] + self.loop[iss].av_corr[t]
  
       for t in range(0,self.tend):
              self.fdis_corr[t] /= self.ns 

       print ("Computed full disconnected correlators for " , self.tag)

   def compute_full_discorr_unbias(self):
       """
          Compute the correlators on each comfiguration and then average
       """

       self.fdis_corr = np.zeros( (self.tend)  )

       for iss in range(0, self.ns) :
         self.loop[iss].compute_corr_unbias(self.tend) 
#         print ("DEBUG corr[0] " , self.loop[iss].av_corr[0] , self.loop[iss].av_corr[1] )
         for t in range(0,self.tend):
              self.fdis_corr[t] = self.fdis_corr[t] + self.loop[iss].av_corr[t]
  
       for t in range(0,self.tend):
              self.fdis_corr[t] /= self.ns 
##              print ("debug " , t , self.fdis_corr[t] ) 

       print ("Computed (UNBIASED) full disconnected correlators for " , self.tag)


   def compute_full_discorr_vacsub(self):

       self.fdis_corr = np.zeros( (self.tend)  )

       vac = np.mean(self.floop)

       for t in range(0,self.tend):
         self.fdis_corr[t] = 0.0
         for tt in range(0,self.nt) :
           tinc = t + tt
           if( tinc >= self.nt):
             tinc = tinc - self.nt
           self.fdis_corr[t] = self.fdis_corr[t] + (self.floop[tt]-vac)*(self.floop[tinc] - vac)
           dd =  (self.floop[t]-vac)*(self.floop[tinc] - vac)

         self.fdis_corr[t] /=  self.nt


   def compute_discorr_err(self):

       self.dis_corr_err  = np.zeros( (self.tend)  )
       tmp = np.zeros( (self.ns)  )
       self.jcorr = np.zeros( (self.nt, self.ns)  )

       for t in range(0,self.tend):
           for i_ss in range(0,self.ns) :
              tmp[i_ss] = 0.0 

              for j_ss in range(0,self.ns) :
                if i_ss != j_ss :
                   tmp[i_ss] += self.loop[j_ss].av_corr[t]

              tmp[i_ss] /=  (self.ns -1 )
              self.jcorr[t][i_ss] = tmp[i_ss]

           self.dis_corr_err[t] = jackknife.jackknife(tmp,self.ns )

       print ("Computed jackknife disconnected correlators for " , self.tag)

   def compute_discorr_vacsub_err(self):

       self.dis_corr_err  = np.zeros( (self.tend)  )
       tmp = np.zeros( (self.ns)  )
       jvac = np.zeros( (self.ns)  )

       for iss in range(0,self.ns) :
          jvac[iss] = 0.0 
          for t in range(0,self.nt):
            jvac[iss] +=  self.jloop[t][iss]
          jvac[iss] /= self.nt


       for t in range(0,self.tend):
           for iss in range(0,self.ns) :
              tmp[iss] = 0.0 
              for tt in range(0,self.nt) :
                tinc = t + tt
                if( tinc >= self.nt):
                   tinc = tinc - self.nt
                tmp[iss] += (self.jloop[tt][iss] -jvac[iss]  ) * (self.jloop[tinc][iss] - jvac[iss] ) 
              tmp[iss] /= self.nt

           self.dis_corr_err[t] = jackknife.jackknife(tmp,self.ns )


   def write_discorr(self):
      """
      Write out the disconnected correlators with errors to standard output
      """

      print ("Disconnected correlators for " , self.tag)
      for t in range(0,self.tend):
          print ("SINFO", t , "%e" % self.fdis_corr[t]   , "%e" % self.dis_corr_err[t] )


   def estimate_amu(self, fname, Vspace):
      """
      Estimate the a_\mu^{HVP (LO) DISC} using the UKQCD/RBC method  (1610.04603).
      """

      print ("Estimate a_mu^{HVP LO}  for " , self.tag)

      print ("Reloading the coefficients from " , fname )
      if not os.path.isfile(fname) :
         print ("Can not find " , fname)
         sys.exit(1)

      print ("SINFO Norm = (taste =4 ) * (spatial volume) * (charge = 9 (3^2))   "  )
      print ("SINFO spatial volume = " , Vspace)
      print ("SINFO T = " , self.tend)

      norm = Vspace * 4 * 9  

      cc = np.load(fname) 
      t = 0 
      for c in cc: 
         print ("w(" , t , ") = " , c )
         t += 1 

      if len(cc) < self.tend :
         print ("Not enough coefficients" )
         sys.exit(1)

      ans = 0.0 
      for t in range(0,self.tend):
         ans +=  self.fdis_corr[t] * cc[t]          

      ans /= norm   


      #  jackknife analysis
      j_amu = np.zeros( (self.ns)  )
      for i_ss in range(0,self.ns) :
         j_amu[i_ss] = 0.0 
         for t in range(0,self.tend):
            for j_ss in range(0,self.ns) :
                if i_ss != j_ss :
                   j_amu[i_ss] += self.loop[j_ss].av_corr[t] * cc[t]

         j_amu[i_ss]  /= (self.ns - 1)

      ans_err = jackknife.jackknife(j_amu,self.ns )
      ans_err /= norm 

      print ("SINFO a_mu^{DIS HVP LO}" ,  "%e" % ans , " +/- %e " % ans_err  )




   def write_dis_loop(self, ttt):
      """
      Write out the disconnected loops to standard output

      """
      print (ttt + "disconnected loop for " , self.tag )
      for t in range(0,self.tend):
          print (t , "%e" % self.floop[t] )


   def save_corr_to_disk(self,fname):
     """
       Save the correlators to disk for each configuration.

     """
     sum = corrutil.corr_sum(self.nt, self.ns, self.tend, self.tag)
     sum.loadin(self.fdis_corr, self.dis_corr_err)

     # write the correlators for each configuration
     for t in range(0,self.tend):
           for iss in range(0,self.ns) :
               sum.all_corr[t][iss] = self.loop[iss].av_corr[t]

     print ("Saving the correlators to the file " , fname )

     f = open(fname, "wb")
     pickle.dump(sum, f)
     f.close()


#
#
#

class dis_corr_3av:
   """
    Class to create disconnected correlators from loop operators.
  There is an average over three components.

   """
   def __init__(self,nt,ns,tend, tag_1, tag_2, tag_3):
       self.ns = 0
       self.nt = nt
       self.tend = tend
       self.tag = tag_1 + "_" + tag_2 + " " + tag_3 
       self.corr_1 = dis_corr(nt,ns,tend, tag_1)
       self.corr_2 = dis_corr(nt,ns,tend, tag_2)
       self.corr_3 = dis_corr(nt,ns,tend, tag_3)



   def load(self, xx_1, xx_2, xx_3):
        """
         Load in the three components of the loop diagram into the class
        """

        self.corr_1.load(xx_1)
        self.corr_2.load(xx_2)
        self.corr_3.load(xx_3)

        self.ns = self.ns + 1

   def  compute_jack_loop(self ):
       """
         Compute the jackknife samples of the loop operators, as
         the average of the three components.

       """
       print ("Computing jackknife loops " , self.tag )
       self.jloop = np.zeros( (self.nt, self.ns)  )

       self.corr_1.compute_jack_loop()
       self.corr_2.compute_jack_loop()
       self.corr_3.compute_jack_loop()

       for t in range(0,self.nt):
         for iss in range(0,self.ns) :
            self.jloop[t][iss] = (self.corr_1.jloop[t][iss] + self.corr_2.jloop[t][iss] + self.corr_3.jloop[t][iss]) / 3.0 



   def  compute_full_loop(self ):

       print ( "Computing full sample loops " , self.tag )

       self.corr_1.compute_full_loop()
       self.corr_2.compute_full_loop()
       self.corr_3.compute_full_loop()


       self.floop = np.zeros( (self.nt)  )
       for t in range(0,self.nt):
         self.floop[t] =  (self.corr_1.floop[t] + self.corr_2.floop[t] + self.corr_3.floop[t] ) / 3.0 


   def dump_full_loop(self ):

       print ( "Full loop " )
       for t in range(0,self.nt):
           print ( "DISCORR_F_DUMP" , t , self.floop[t] )
   
   def compute(self):
       """
          Compute the full and jacksamples of the loops.
          Use the loops to compute the disconnected correlators.

       """

       self.compute_full_loop()
       self.compute_jack_loop()

       self.compute_full_discorr() 
       self.compute_discorr_err()

   def compute_unbias(self):
       """
          Compute the full and jacksamples of the loops.
          Use the loops to compute the disconnected correlators.

       """

       self.compute_full_loop()
       self.compute_jack_loop()

       self.compute_full_discorr_unbias() 
       self.compute_discorr_err()


   def compute_full_discorr(self):
       """
          Compute the full sample disconnected correlators

       """


       self.corr_1.compute_full_discorr()
       self.corr_2.compute_full_discorr()
       self.corr_3.compute_full_discorr()

       self.fdis_corr = np.zeros( (self.tend)  )

       for t in range(0,self.tend):
         self.fdis_corr[t] = ( self.corr_1.fdis_corr[t] + self.corr_2.fdis_corr[t] + self.corr_3.fdis_corr[t] ) / 3.0 



   def compute_full_discorr_unbias(self):
       """
          Compute the full sample disconnected correlators

       """

       self.corr_1.compute_full_discorr_unbias()
       self.corr_2.compute_full_discorr_unbias()
       self.corr_3.compute_full_discorr_unbias()

       self.fdis_corr = np.zeros( (self.tend)  )

       for t in range(0,self.tend):
         self.fdis_corr[t] = ( self.corr_1.fdis_corr[t] + self.corr_2.fdis_corr[t] + self.corr_3.fdis_corr[t] ) / 3.0 



   def compute_discorr_err(self):

       self.dis_corr_err  = np.zeros( (self.tend)  )
       tmp = np.zeros( (self.ns)  )

       self.corr_1.compute_discorr_err()
       self.corr_2.compute_discorr_err()
       self.corr_3.compute_discorr_err()

       for t in range(0,self.tend):
           for i_ss in range(0,self.ns) :
              tmp[i_ss] = (self.corr_1.jcorr[t][i_ss] + self.corr_2.jcorr[t][i_ss] + self.corr_3.jcorr[t][i_ss]) /3.0 

           self.dis_corr_err[t] = jackknife.jackknife(tmp,self.ns )


   def write_discorr(self):
      """
      Write out the disconnected correlators with errors to standard output
      """

      print ("Disconnected correlators for " , self.tag )
      for t in range(0,self.tend):
          print ("SINFO", t , "%e" % self.fdis_corr[t]   , "%e" % self.dis_corr_err[t] )


   def write_dis_loop(self):
      """
      Write out the disconnected loops to standard output

      """
      print ("Disconnected loop for " , self.tag )
      for t in range(0,self.tend):
          print (t , "%e" % self.floop[t] )


   def save_corr_to_disk(self,fname):
     """
     Save the correlators to disk for each configuration. 
     """

     sum = corrutil.corr_sum(self.nt, self.ns, self.tend, self.tag)
     sum.loadin(self.fdis_corr, self.dis_corr_err)

     # write the correlators for each configuration
     for t in range(0,self.tend):
           for iss in range(0,self.ns) :
               sum.all_corr[t][iss] = (self.corr_1.loop[iss].av_corr[t] + self.corr_2.loop[iss].av_corr[t] + self.corr_3.loop[iss].av_corr[t] ) / 3.0 

     print ("Saving the correlators to the file " , fname )

     f = open(fname, "wb")
     pickle.dump(sum, f)
     f.close()



   def estimate_amu(self, fname, Vspace):
       print ( "estimate_amu::Not implemented" )



####

class dis_corr_isospin:
   """
   Class to create disconnected correlators from loop operators with an estimate
   of isospin breaking.

   """

   def __init__(self,nt,ns,tend, tag):
       self.ns = 0
       self.nt = nt
       self.tend = tend
       self.tag = tag 
       self.loop_us = [] 
       self.loop_ud = [] 


   def load(self, xx_us, xx_ud):
        """
          Load in the loops 
        """
        self.loop_us.append(xx_us)
        self.loop_ud.append(xx_ud)

        self.ns = self.ns + 1

   
   def compute(self, mass_us):
       """
          Use the loops to compute the disconnected correlators.

       """
       self.compute_full_discorr(mass_us) 
       self.compute_discorr_err()

   def compute_unbias(self, mass_us):
       """
          Use the loops to compute the disconnected correlators.
          (Remove the single noise sources)

       """
       self.compute_full_discorr_unbias(mass_us) 
       self.compute_discorr_err()


   def compute_off(self, mass_us):
       """
          Use the loops to compute the disconnected correlators.
          Off diagonal comtribution
       """
       self.compute_full_discorr_off(mass_us) 
       self.compute_discorr_err()


   def compute_full_discorr(self,  mass_ud):
       """
          Compute the correlators on each comfiguration and then average
       """

       self.fdis_corr = np.zeros( (self.tend)  )

       for iss in range(0, self.ns) :
         self.loop_us[iss].compute_corr_isospin(self.tend, self.loop_ud[iss], mass_ud) 
         for t in range(0,self.tend):
              self.fdis_corr[t] = self.fdis_corr[t] + self.loop_us[iss].av_corr[t]
  
       for t in range(0,self.tend):
              self.fdis_corr[t] /= self.ns 

       print ( "Computed full ISOSPIN disconnected correlators for " , self.tag, " delta_mass = " , mass_ud )



   def compute_full_discorr_unbias(self,  mass_ud):
       """
          Compute the correlators on each comfiguration and then average
       """

       self.fdis_corr = np.zeros( (self.tend)  )

       for iss in range(0, self.ns) :
         self.loop_us[iss].compute_corr_isospin_unbias(self.tend, self.loop_ud[iss], mass_ud) 
         for t in range(0,self.tend):
              self.fdis_corr[t] = self.fdis_corr[t] + self.loop_us[iss].av_corr[t]
  
       for t in range(0,self.tend):
              self.fdis_corr[t] /= self.ns 

       print ( "Computed full (unbias) ISOSPIN disconnected correlators for " , self.tag, " delta_mass = " , mass_ud )



   def compute_full_discorr_off(self,  mass_ud):
       """
          Compute the correlators on each comfiguration and then average
          Off diagonal contribution.
       """

       self.fdis_corr = np.zeros( (self.tend)  )

       for iss in range(0, self.ns) :
         self.loop_us[iss].compute_corr_isospin_off(self.tend, self.loop_ud[iss], mass_ud) 
         for t in range(0,self.tend):
              self.fdis_corr[t] = self.fdis_corr[t] + self.loop_us[iss].av_corr[t]
  
       for t in range(0,self.tend):
              self.fdis_corr[t] /= self.ns 

       print ( "Computed full ISOSPIN disconnected correlators for " , self.tag, " delta_mass = " , mass_ud )



   def compute_discorr_err(self):

       self.dis_corr_err  = np.zeros( (self.tend)  )
       tmp = np.zeros( (self.ns)  )
       self.jcorr = np.zeros( (self.nt, self.ns)  )

       for t in range(0,self.tend):
           for i_ss in range(0,self.ns) :
              tmp[i_ss] = 0.0 

              for j_ss in range(0,self.ns) :
                if i_ss != j_ss :
                   tmp[i_ss] += self.loop_us[j_ss].av_corr[t]

              tmp[i_ss] /=  (self.ns -1 )
              self.jcorr[t][i_ss] = tmp[i_ss]

           self.dis_corr_err[t] = jackknife.jackknife(tmp,self.ns )

       print ("Computed jackknife isospin disconnected correlators for " , self.tag )


   def write_discorr(self):
      """
      Write out the disconnected correlators with errors to standard output
      """

      print ("Disconnected correlators for " , self.tag )
      for t in range(0,self.tend):
          print ("SINFO isospin", t , "%e" % self.fdis_corr[t]   , "%e" % self.dis_corr_err[t] )


   def estimate_amu(self, fname, Vspace):
      """
      Estimate the a_\mu^{HVP (LO) DISC} using the UKQCD/RBC method  (1610.04603).
      """

      print ("Estimate a_mu^{HVP LO}  for " , self.tag )

      print ("Reloading the coefficients from " , fname  )
      if not os.path.isfile(fname) :
         print ("Can not find " , fname )
         sys.exit(1)

      print ("SINFO Norm = (taste =4 ) * (spatial volume) * (charge = 1/9)   " , )
      print ("SINFO spatial volume = " , Vspace )
      print ("SINFO T = " , self.tend)

      norm = Vspace * 4 * 9  

      cc = np.load(fname) 
      t = 0 
      for c in cc: 
         print ("w(" , t , ") = " , c )
         t += 1 

      if len(cc) < self.tend :
         print ("Not enough coefficients" )
         sys.exit(1)

      ans = 0.0 
      for t in range(0,self.tend):
         ans +=  self.fdis_corr[t] * cc[t]          

      ans /= norm   


      #  jackknife analysis
      j_amu = np.zeros( (self.ns)  )
      for i_ss in range(0,self.ns) :
         j_amu[i_ss] = 0.0 
         for t in range(0,self.tend):
            for j_ss in range(0,self.ns) :
                if i_ss != j_ss :
                   j_amu[i_ss] += self.loop_us[j_ss].av_corr[t] * cc[t]

         j_amu[i_ss]  /= (self.ns - 1)

      ans_err = jackknife.jackknife(j_amu,self.ns )
      ans_err /= norm 

      print ("SINFO a_mu^{DIS HVP LO}" ,  "%e" % ans , " +/- %e " % ans_err  )




   def save_corr_to_disk(self,fname):
     """
       Save the correlators to disk for each configuration.

     """
     sum = corrutil.corr_sum(self.nt, self.ns, self.tend, self.tag)
     sum.loadin(self.fdis_corr, self.dis_corr_err)

     # write the correlators for each configuration
     for t in range(0,self.tend):
           for iss in range(0,self.ns) :
               sum.all_corr[t][iss] = self.loop_us[iss].av_corr[t]

     print ("Saving the correlators to the file " , fname )

     f = open(fname, "wb")
     pickle.dump(sum, f)
     f.close()





####

class dis_corr_isospin_3av:
   """
   Class to create disconnected correlators from loop operators with an estimate
   of isospin breaking.  Use all three components.

   """

   def __init__(self,nt,ns,tend, tag_1):
       self.ns = 0
       self.nt = nt
       self.tend = tend
       self.tag_1 = tag_1 
       self.loop_1_us = [] 
       self.loop_1_ud = [] 

       self.loop_2_us = [] 
       self.loop_2_ud = [] 

       self.loop_3_us = [] 
       self.loop_3_ud = [] 


   def load(self, xx_1_us, xx_1_ud, xx_2_us, xx_2_ud, xx_3_us, xx_3_ud):
        """
          Load in the loops 
        """
        self.loop_1_us.append(xx_1_us)
        self.loop_1_ud.append(xx_1_ud)

        self.loop_2_us.append(xx_2_us)
        self.loop_2_ud.append(xx_2_ud)

        self.loop_3_us.append(xx_3_us)
        self.loop_3_ud.append(xx_3_ud)

        self.ns = self.ns + 1

   
   def compute(self, mass_us):
       """
          Use the loops to compute the disconnected correlators.

       """
       self.compute_full_discorr(mass_us) 
       self.compute_discorr_err()

   def compute_unbias(self, mass_us):
       """
          Use the loops to compute the disconnected correlators.
          (Remove the single noise sources)

       """
       self.compute_full_discorr_unbias(mass_us) 
       self.compute_discorr_err()


   def compute_off(self, mass_us):
       """
          Use the loops to compute the disconnected correlators.
          Off diagonal comtribution
       """
       self.compute_full_discorr_off(mass_us) 
       self.compute_discorr_err()


   def compute_full_discorr(self,  mass_ud):
       """
          Compute the correlators on each comfiguration and then average
       """

       self.fdis_corr = np.zeros( (self.tend)  )

       for iss in range(0, self.ns) :
         self.loop_1_us[iss].compute_corr_isospin(self.tend, self.loop_1_ud[iss], mass_ud) 
         self.loop_2_us[iss].compute_corr_isospin(self.tend, self.loop_2_ud[iss], mass_ud) 
         self.loop_3_us[iss].compute_corr_isospin(self.tend, self.loop_3_ud[iss], mass_ud) 

         for t in range(0,self.tend):
              self.fdis_corr[t] = self.fdis_corr[t] + (self.loop_1_us[iss].av_corr[t] + self.loop_2_us[iss].av_corr[t] + self.loop_3_us[iss].av_corr[t]) / 3.0
  
       for t in range(0,self.tend):
              self.fdis_corr[t] /= self.ns 

       print ( "Computed full ISOSPIN disconnected correlators for " , self.tag_1, " delta_mass = " , mass_ud )



   def compute_full_discorr_unbias(self,  mass_ud):
       """
          Compute the correlators on each comfiguration and then average
       """

       self.fdis_corr = np.zeros( (self.tend)  )

       for iss in range(0, self.ns) :
         self.loop_1_us[iss].compute_corr_isospin_unbias(self.tend, self.loop_1_ud[iss], mass_ud) 
         self.loop_2_us[iss].compute_corr_isospin_unbias(self.tend, self.loop_2_ud[iss], mass_ud) 
         self.loop_3_us[iss].compute_corr_isospin_unbias(self.tend, self.loop_3_ud[iss], mass_ud) 

         for t in range(0,self.tend):
              self.fdis_corr[t] = self.fdis_corr[t] + (self.loop_1_us[iss].av_corr[t] + self.loop_2_us[iss].av_corr[t] + self.loop_3_us[iss].av_corr[t]) / 3.0 
  
       for t in range(0,self.tend):
              self.fdis_corr[t] /= self.ns 

       print ( "Computed full (unbias) ISOSPIN disconnected correlators for " , self.tag_1, " delta_mass = " , mass_ud )



   def compute_full_discorr_off(self,  mass_ud):
       """
          Compute the correlators on each comfiguration and then average
          Off diagonal contribution.
       """

       self.fdis_corr = np.zeros( (self.tend)  )

       for iss in range(0, self.ns) :
         self.loop_1_us[iss].compute_corr_isospin_off(self.tend, self.loop_1_ud[iss], mass_ud) 
         self.loop_2_us[iss].compute_corr_isospin_off(self.tend, self.loop_2_ud[iss], mass_ud) 
         self.loop_3_us[iss].compute_corr_isospin_off(self.tend, self.loop_3_ud[iss], mass_ud) 

         for t in range(0,self.tend):
              self.fdis_corr[t] = self.fdis_corr[t] + ( self.loop_1_us[iss].av_corr[t] +  self.loop_2_us[iss].av_corr[t] +  self.loop_3_us[iss].av_corr[t] ) / 3.0
  
       for t in range(0,self.tend):
              self.fdis_corr[t] /= self.ns 

       print ( "Computed full ISOSPIN disconnected correlators for " , self.tag_1, " delta_mass = " , mass_ud )



   def compute_discorr_err(self):

       self.dis_corr_err  = np.zeros( (self.tend)  )
       tmp = np.zeros( (self.ns)  )
       self.jcorr = np.zeros( (self.nt, self.ns)  )

       for t in range(0,self.tend):
           for i_ss in range(0,self.ns) :
              tmp[i_ss] = 0.0 

              for j_ss in range(0,self.ns) :
                if i_ss != j_ss :
                   tmp[i_ss] += (self.loop_1_us[j_ss].av_corr[t] + self.loop_2_us[j_ss].av_corr[t] + self.loop_3_us[j_ss].av_corr[t]) / 3.0

              tmp[i_ss] /=  (self.ns -1 )
              self.jcorr[t][i_ss] = tmp[i_ss]

           self.dis_corr_err[t] = jackknife.jackknife(tmp,self.ns )

       print ("Computed jackknife isospin disconnected correlators for " , self.tag_1 )


   def write_discorr(self):
      """
      Write out the disconnected correlators with errors to standard output
      """

      print ("Disconnected correlators for " , self.tag_1 )
      for t in range(0,self.tend):
          print ("SINFO isospin", t , "%e" % self.fdis_corr[t]   , "%e" % self.dis_corr_err[t] )


   def estimate_amu(self, fname, Vspace):
      """
      Estimate the a_\mu^{HVP (LO) DISC} using the UKQCD/RBC method  (1610.04603).
      """

      print ("Estimate a_mu^{HVP LO}  for " , self.tag_1 )

      print ("Reloading the coefficients from " , fname  )
      if not os.path.isfile(fname) :
         print ("Can not find " , fname )
         sys.exit(1)

      print ("SINFO Norm = (taste =4 ) * (spatial volume) * (charge = 1/9)   " , )
      print ("SINFO spatial volume = " , Vspace )
      print ("SINFO T = " , self.tend)

      norm = Vspace * 4 * 9  

      cc = np.load(fname) 
      t = 0 
      for c in cc: 
         print ("w(" , t , ") = " , c )
         t += 1 

      if len(cc) < self.tend :
         print ("Not enough coefficients" )
         sys.exit(1)

      ans = 0.0 
      for t in range(0,self.tend):
         ans +=  self.fdis_corr[t] * cc[t]          

      ans /= norm   


      #  jackknife analysis
      j_amu = np.zeros( (self.ns)  )
      for i_ss in range(0,self.ns) :
         j_amu[i_ss] = 0.0 
         for t in range(0,self.tend):
            for j_ss in range(0,self.ns) :
                if i_ss != j_ss :
                   j_amu[i_ss] += (self.loop_1_us[j_ss].av_corr[t] + self.loop_2_us[j_ss].av_corr[t] + self.loop_3_us[j_ss].av_corr[t] ) / 3.0     * cc[t]

         j_amu[i_ss]  /= (self.ns - 1)

      ans_err = jackknife.jackknife(j_amu,self.ns )
      ans_err /= norm 

      print ("SINFO a_mu^{DIS HVP LO}" ,  "%e" % ans , " +/- %e " % ans_err  )




   def save_corr_to_disk(self,fname):
     """
       Save the correlators to disk for each configuration.

     """
     sum = corrutil.corr_sum(self.nt, self.ns, self.tend, self.tag_1)
     sum.loadin(self.fdis_corr, self.dis_corr_err)

     # write the correlators for each configuration
     for t in range(0,self.tend):
           for iss in range(0,self.ns) :
               sum.all_corr[t][iss] = (self.loop_1_us[iss].av_corr[t] + self.loop_2_us[iss].av_corr[t] + self.loop_3_us[iss].av_corr[t] ) / 3.0

     print ("Saving the correlators to the file " , fname )

     f = open(fname, "wb")
     pickle.dump(sum, f)
     f.close()


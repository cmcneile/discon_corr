import re
import math
import numpy as np
import sys


""""
    Class for reading the stochastic estimates of disconnected loops
"""



class single_loop:
   """
      Class to store estimates of the quark loops on a single configuration,
      Stochastic estimate of the loops.
   """
   def __init__(self,nt,ns,tag):
       self.nt = nt 
       self.ns = ns
       self.STATns = ns
       self.tag = tag
       self.mass = 0
       self.loop = np.zeros( (nt, ns)  ) 
#       self.loop_tmp  = np.zeros( ( ns)  ) 
#       self.loop_tmpA  = np.zeros( ( nt)  ) 
       self.av_loop    = np.zeros( (nt)  ) 
       self.av_corr    = np.zeros( (nt)  ) 
       self.av_loop_err = np.zeros( (nt)  ) 
       self.noise_pt = 8

##       self.loop_pt = 5
       self.loop_pt = 6

       if re.search(tag , "DIS_LOOP_one_hpqcd") :
         self.noise_pt = 7
         self.loop_pt = 5


       print ("INFO Noise pointer = " , self.noise_pt)
       print ("INFO Correlator pointer = " , self.loop_pt)


   def set_re(self):
      print ("Using real loop")
      self.loop_pt = 5


   def load_data(self,filename):
      tag = self.tag
      print ("Reading from " , filename , " searching for " , tag) 

      f = open(filename, 'rU')
      text = f.readlines()

      count = 0 
      
      for line in text:
         lll = line.rstrip()
         if re.search(tag,  lll):
#            print lll
#            print tmp[4] , tmp[8] , tmp[6]
            tmp = lll.split()
            ccc = float(tmp[self.loop_pt])


            ii = int(tmp[ self.noise_pt]) - 1
            tt = int(tmp[4])
            self.loop[tt][ii] = ccc

      f.close()


   def load_data_mass(self,filename, mass_in):
      tag = self.tag
      print ("Reading from " , filename , " searching for " , tag, " mass = " , mass_in )
      self.mass = mass_in

      f = open(filename, 'rU')
      text = f.readlines()

      count = 0 
      
      for line in text:
         lll = line.rstrip()
         if re.search(tag,  lll):
#            print lll
#            print tmp[4] , tmp[8] , tmp[6]
            tmp = lll.split()
            mass = float(tmp[2] )
            if mass == mass_in :  
              ccc = float(tmp[self.loop_pt])
              ii = int(tmp[ self.noise_pt]) - 1
              tt = int(tmp[4])
              self.loop[tt][ii] = ccc

      f.close()



   def load_data_mass_more(self,filename, mass_in, noise_start):
      tag = self.tag
      print ("Reading from " , filename , " searching for " , tag, " mass = " , mass_in )
      self.mass = mass_in

      f = open(filename, 'rU')
      text = f.readlines()

      count = 0 
      
      for line in text:
         lll = line.rstrip()
         if re.search(tag,  lll):
#            print lll
#            print tmp[4] , tmp[8] , tmp[6]
            tmp = lll.split()
            mass = float(tmp[2] )
            if mass == mass_in :  
              ccc = float(tmp[self.loop_pt])
              ii = int(tmp[ self.noise_pt]) - 1  + noise_start
              tt = int(tmp[4])
              self.loop[tt][ii] = ccc

      f.close()


   def load_data_inc(self,filename,tag):
      print ("Reading (inc) from " , filename , " searching for " , tag )

      f = open(filename, 'rU')
      text = f.readlines()

      count = 0 
      
      for line in text:
         lll = line.rstrip()
         if re.search(tag,  lll):
#            print lll
#            print tmp[4] , tmp[8] , tmp[6]
            tmp = lll.split()
            ccc = float(tmp[self.loop_pt])


            ii = int(tmp[ self.noise_pt]) - 1
            tt = int(tmp[4])
            self.loop[tt][ii] += ccc

      f.close()



   def load_data_mass_inc(self,filename,tag, mass_in):
      print ("Reading (inc) from " , filename , " searching for " , tag, "mass= " , mass_in )

      f = open(filename, 'rU')
      text = f.readlines()

      count = 0 
      
      for line in text:
         lll = line.rstrip()
         if re.search(tag,  lll):
#            print lll
#            print tmp[4] , tmp[8] , tmp[6]
            tmp = lll.split()
            mass = float(tmp[2] )
            if mass == mass_in :  
              ccc = float(tmp[self.loop_pt])
              ii = int(tmp[ self.noise_pt]) - 1
              tt = int(tmp[4])
              self.loop[tt][ii] += ccc

      f.close()





   def load_data_inc_old(self,filename,tag):
      print ("Reading from " , filename , " searching for " , tag )

      f = open(filename, 'rU')
      text = f.readlines()


      count = 0 
      
      for line in text:
         lll = line.rstrip()
         if re.search(tag,  lll):
#            print lll
            tmp = lll.split()
            ccc = float(tmp[1])
#            print ccc

            dd = tmp[0] 
            ddd = dd.replace(tag,"")
            dd  = ddd.replace("]["," ")
            ddd = dd.replace("]","")
            dd  = ddd.replace("[","")
            eee = dd.split()

            ii = int(eee[0]) -1
            tt = int(eee[1])
#            print int(eee[0]) , int(eee[1])
#            print tt , int(eee[0]) 
            self.loop[tt][ii] = self.loop[tt][ii] + ccc

      f.close()


##
##

   def dump_time_avloop(self):

      print ("Dump for " , self.tag )
      for tt in range(0,self.nt):
         print ( tt , self.av_loop[tt], self.av_loop_err[tt] )

   def dump_loop(self):

      print ( "Dump loop for " , self.tag , "mass = " , self.mass )
      for tt in range(0,self.nt):
         for ii in range(0,self.STATns):
              print ( tt , ii, "%e " % self.loop[tt][ii] )
  

##
##

   def norm_corr(self, norm):

      print ("Normalising correction data with factor 1/ " , norm )
      for tt in range(0,self.nt):
        self.av_loop[tt]     = self.av_loop[tt] / norm
        self.av_loop_err[tt] = self.av_loop_err[tt] / norm

      self.loop_t_err = self.loop_t_err / norm


##
##

   def dump_summary(self,ccc,norm=1):

      print ("Summary for " , self.tag )
      if norm != 1 :
         print ("Normalisation eigen " , norm )

      print ("Time\t Eigen\t\tcorrection\t\t\tTotal\t\t err" )
      for tt in range(0,self.nt):
         ccc_l = ccc[tt] / norm 
         total = ccc_l + self.av_loop[tt]
         print (tt , "\t" ,  )
         print ("%8.5f" , ccc_l , )
         print ("\t" , self.av_loop[tt], "\t" , total , "\t" , self.av_loop_err[tt] )


##
##

   def dump_summary_cmp(self,ccc,localCorr,norm=1):

      print ("Summary for " , self.tag )
      if norm != 1 :
         print ("Normalisation eigen " , norm )

      print ("Time\t Eigen + \tcorrection = \tTotal\t\t\tLocal_point_corr" )
      for tt in range(0,self.nt):
         ccc_l = ccc[tt] / norm 
         total = ccc_l + self.av_loop[tt]
         print ("CORR[", tt , "]\t" ,  )
         print ("%8.6f" % ccc_l ,  )
         print ("\t" , "%8.6f" %  self.av_loop[tt], )
         print ("\t" , "%8.6f" % total  , "  +/- " , )
         print   (     "%8.6f" % self.av_loop_err[tt],  )
         print (" lpion " , "%8.6f" %  localCorr[tt] ,  )
         sig = ( localCorr[tt] - total ) / self.av_loop_err[tt]
         print (" no. sig " , sig )




##
##

   def add_corr(self,ccc,norm=1):

      print ("Summary for " , self.tag )
      if norm != 1 :
         print ("add_corr:  Normalisation eigen " , norm )

      for tt in range(0,self.nt):
         self.av_loop[tt] = self.av_loop[tt] + ccc[tt] / norm 


   def write_to_disk(self,fname) :

     file = open(fname,"w")
     for tt in range(0,self.nt):
       print >> file , ("%e" %  self.av_loop[tt] )

##       print >> file , tt, self.av_loop[tt]

     file.close()



##
##

   def dump_summary_cmpA(self,ccc,localCorr,localCorr_err,norm=1):

      print ("Summary for " , self.tag )
      if norm != 1 :
         print ("Normalisation eigen " , norm )

      print ("Time\t Eigen + \tcorrection = \tTotal\t\t\tLocal_point_corr" )
      for tt in range(0,self.nt):
         ccc_l = ccc[tt] / norm 
         total = ccc_l + self.av_loop[tt]
         print ( "CORR[", tt , "]\t" ,  )
         print ( "%8.6f" % ccc_l ,  )
         print ( "\t" , "%8.6f" %  self.av_loop[tt],  )
         print ( "\t" , "%8.6f" % total  , "  +/- " , )
         print      (  "%8.6f" % self.av_loop_err[tt], )
         print ( " lpion " , "%8.6f" %  localCorr[tt] ,  )
         sig = ( localCorr[tt] - total ) / localCorr_err[tt]
         print ( " no. sig " , sig )



   def dump_summary_tsum_cmp(self,ccc,localCorr,ne,norm=1):

      print ( "Summary for " , self.tag )
      if norm != 1 :
         print ( "Normalisation eigen " , norm )

      s_ccc_l   = 0.0 
      s_total   = 0.0
      s_av_loop = 0.0
      s_pt      = 0.0

      for tt in range(0,self.nt):
         ccc_l = ccc[tt] / norm 
         s_ccc_l = s_ccc_l + ccc_l

         total = ccc_l + self.av_loop[tt]
         s_total = s_total + total
         s_av_loop = s_av_loop + self.av_loop[tt]

         s_pt = s_pt + localCorr[tt]


##      print "Eigen + \tcorrection = \t\tTotal\t\t\tLocal_point_corr"
      print ( "Eigen + \tcorrection = \tTotal" )
##      print , 
      print (ne, )
      print ( "CORR %8.6f" % s_ccc_l ,  )
      print ( "\t" , "%8.6f" %  s_av_loop,  )

      print   ( "\t%8.6f" %  s_total,  )
      tmp_err = self.loop_t_err 
      print   ( "\t%8.6f" % tmp_err , )

      print   ( "\t%8.6f" %  s_pt )

##      print "DEBUG " , self.loop_t , " +/- " , self.loop_t_err


   def adjust_samples_compute(self,ns) :
      self.STATns = ns

##
##

   def compute(self):
      """"
       Average over the noise sources to estimate the loops on each time slice
      """

      corr_tmp = np.zeros( ( self.STATns  ) )
      print ( "Compute average correlators for " , self.tag, " using nsample" , self.STATns )
      for tt in range(0,self.nt):
         corr_re = 0.0 
         for ii in range(0,self.STATns):
            corr_tmp[ii] = self.loop[tt][ii]

         corr_re = np.mean(corr_tmp)
         corr_err = np.std(corr_tmp) / math.sqrt(self.STATns) 
         self.av_loop[tt] = corr_re
         self.av_loop_err[tt] = corr_err
#         print tt , corr_re, corr_err


   def compute_corr(self,tend) :
      """
         Compute the correlators from the loops
         This estimate includes the bias from equal noise sources.
      """

      for t in range(0,self.nt):
         self.av_corr[t] = 0.0
         for tt in range(0,self.nt) :
           tinc = t + tt
           if( tinc >= self.nt):
             tinc = tinc - self.nt
           self.av_corr[t] = self.av_corr[t] + self.av_loop[tt]*self.av_loop[tinc]
         self.av_corr[t] /=  self.nt


   def compute_corr_unbias(self,tend) :
      """
         Compute the correlators from the loops
         This estimate removes the equal time bias
      """

      for t in range(0,self.nt):
         self.av_corr[t] = 0.0
         for tt in range(0,self.nt) :
           tinc = t + tt
           if( tinc >= self.nt):
             tinc = tinc - self.nt
           self.av_corr[t] = self.av_corr[t] + self.av_loop[tt]*self.av_loop[tinc]
           bias = 0.0 
           for ii in range(0,self.STATns):
              bias += self.loop[tt][ii] * self.loop[tinc][ii]
##           print ("BIAS " , t , tt , bias , " corr= " , self.av_corr[t]  )   
           self.av_corr[t] -= bias /(self.STATns**2)

         self.av_corr[t] /=  self.nt


   def compute_corr_isospin(self,tend, loop_ud, mass_ud) :
      """
         Compute the isospin correlators from the loops
         This includes the contamination from equal noise sources
      """

      for t in range(0,self.nt):
         self.av_corr[t] = 0.0
         for tt in range(0,self.nt) :
           tinc = t + tt
           if( tinc >= self.nt):
             tinc = tinc - self.nt
           tmp = (self.av_loop[tt] + 2.0*mass_ud*loop_ud.av_loop[tt] ) * ( self.av_loop[tinc] + 2.0*mass_ud * loop_ud.av_loop[tinc] )
           self.av_corr[t] = self.av_corr[t] + tmp  

         self.av_corr[t] /=  self.nt


   def compute_corr_isospin_unbias(self,tend, loop_ud, mass_ud) :
      """
         Compute the isospin correlators from the loops
         This includes the contamination from equal noise sources
      """

      for t in range(0,self.nt):
         self.av_corr[t] = 0.0
         for tt in range(0,self.nt) :
           tinc = t + tt
           if( tinc >= self.nt):
             tinc = tinc - self.nt
           tmp = (self.av_loop[tt] + 2.0 * mass_ud * loop_ud.av_loop[tt] ) * ( self.av_loop[tinc] + 2.0 * mass_ud * loop_ud.av_loop[tinc] )
           self.av_corr[t] = self.av_corr[t] + tmp  
           bias = 0.0 
           for ii in range(0,self.STATns):
              bias += (self.loop[tt][ii] + 2.0*mass_ud*loop_ud.loop[tt][ii] ) * (self.loop[tinc][ii]  + 2.0*mass_ud * loop_ud.loop[tinc][ii]  )


##           print ("BIAS " , t , tt , bias , " corr= " , self.av_corr[t]  )   
           self.av_corr[t] -= bias /(self.STATns**2)


         self.av_corr[t] /=  self.nt



   def compute_corr_isospin_off(self,tend, loop_ud, mass_ud) :
      """
         Compute the isospin correlators from the loops
         This is the off diagonl contribution
      """

      for t in range(0,self.nt):
         self.av_corr[t] = 0.0
         for tt in range(0,self.nt) :
           tinc = t + tt
           if( tinc >= self.nt):
             tinc = tinc - self.nt
           tmp = (self.av_loop[tt]  ) * ( 2.0 *  mass_ud * loop_ud.av_loop[tinc] )
           self.av_corr[t] = self.av_corr[t] + tmp  

         self.av_corr[t] /=  self.nt



   def include_eigen_corr(self,eigen):
      """"
      Add in the loops estimated from the eigenvalues
      """

      print ( "The eigenpair loops are being added" )
      for tt in range(0,self.nt):
         self.av_loop[tt] += eigen.corr[tt]

   def include_eigen_corr_fact(self,eigen,fact):
      """"
      Add in the loops estimated from the eigenvalues
      """

      print ("The eigenpair loops are being added with norm = " , fact )
      for tt in range(0,self.nt):
         self.av_loop[tt] += fact * eigen.corr[tt]



#
#  Replace the correlators with those estimated from the eigenvalues
#
   def replace_with_eigen_corr(self,eigen):

      print ( "The eigenpair correlators REPLACING stochastic correlators" )
      for tt in range(0,self.nt):
#         print "DEBUG " , self.av_loop[tt] ,  eigen.corr[tt]
         self.av_loop[tt] = eigen.corr[tt]


#
#  Add in the correlators estimated from the eigenvalues
#
   def include_eigen_corr_strange(self,strange, eigen, eigen_strange):

      print ( "The eigenpair correlators are being added" )
      for tt in range(0,self.nt):
#         print "DEBUG " , self.av_loop[tt] ,  eigen.corr[tt]


          self.av_loop[tt] += eigen.corr[tt] - (eigen_strange.corr[tt] + strange.av_loop[tt])



#
#  Add in the correlators estimated from the eigenvalues
#
   def subtract(self,strange):
      """
         Subtract the strange loops from the light loops
         Work at the average loops and the lower level loops
      """
      print ( "The light-strange correlators are being created" )
      for tt in range(0,self.nt):
          self.av_loop[tt] -=  strange.av_loop[tt]
          for ii in range(0,self.STATns):
             self.loop[tt][ii] -= strange.loop[tt][ii]

   def include_eigen_corr_strange_A(self,strange, eigen, eigen_strange):

      print ( "The eigenpair correlators are being added" )
      for tt in range(0,self.nt):
#         print "DEBUG " , self.av_loop[tt] ,  eigen.corr[tt]
          self.av_loop[tt] = eigen.corr[tt] - eigen_strange.corr[tt] 

   def include_eigen_corr_strange_B(self,strange, eigen, eigen_strange):

      print ("The eigenpair correlators are being added" )
      for tt in range(0,self.nt):
#         print "DEBUG " , self.av_loop[tt] ,  eigen.corr[tt]


          self.av_loop[tt] -=  strange.av_loop[tt]



#
#  Compute the average of the corr summed over time
#
   def compute_time_sum(self):

     corr_tmp = np.zeros( ( self.STATns  ) )
     for ii in range(0,self.STATns):
       tot = 0 
       for tt in range(0,self.nt):
          tot = tot + self.loop[tt][ii]
       corr_tmp[ii] = tot


     tot = 0 
     for tt in range(0,self.nt):
          tot = tot + self.av_loop[tt]

     self.loop_t = tot

##     self.loop_t     = np.mean(corr_tmp)
     self.loop_t_err = np.std( corr_tmp) / math.sqrt(self.STATns) 




   def combine(self,a,b,c):
      print ( "Combination of 3 components  " , self.tag )
      for tt in range(0,self.nt):
         self.av_loop[tt] = a.av_loop[tt] + b.av_loop[tt] + c.av_loop[tt] 
         self.av_loop_err[tt] = math.sqrt(a.av_loop_err[tt]**2 + b.av_loop_err[tt]**2 + c.av_loop_err[tt]**2 ) 

      self.loop_t = a.loop_t + b.loop_t + c.loop_t
      self.loop_t_err = math.sqrt(a.loop_t_err**2 + b.loop_t_err**2 + c.loop_t_err**2)


   def combine_4(self,a,b,c,d):
      print ( "Combination of 3 components  " , self.tag )
      for tt in range(0,self.nt):
         self.av_loop[tt] = a.av_loop[tt] + b.av_loop[tt] + c.av_loop[tt]  + d.av_loop[tt] 
         self.av_loop_err[tt] = math.sqrt(a.av_loop_err[tt]**2 + b.av_loop_err[tt]**2 + c.av_loop_err[tt]**2 + d.av_loop_err[tt]**2 ) 

      self.loop_t = a.loop_t + b.loop_t + c.loop_t + d.loop_t
      self.loop_t_err = math.sqrt(a.loop_t_err**2 + b.loop_t_err**2 + c.loop_t_err**2 + d.loop_t_err**2)


   def norm_lower(self,norm):

      print ( "Applying norm = " , norm, " to correlators" )
      for tt in range(0,self.nt):
         for ii in range(0,self.ns):
              self.loop[tt][ii] /= norm

   def cmp_time_sum(self,val):
         nosigma = (self.loop_t - val) / self.loop_t_err 
         print ( "INFO TSUM_SIGMA_DEV", self.tag ,  "nosigma= %f" % nosigma , "deviation from %f" % val  )


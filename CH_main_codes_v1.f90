!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccc
!ccc   FORTRAN Codes for time-varying coeff. models & testing the equivalence of hazards   
!ccc   (Bayesian Analysis) 
!ccc
!ccc   2014.01.01 / 2016.04.20 (rearranged) 
!ccc   by G. Kim 
!ccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       use                random   

      integer             od, d, sh, itr, burn, num, loc
      integer, dimension(12) :: seed = (/1,3,5,7,8,10,3,4,5,6,7,8/) 

      double precision    sm, sde1, sde2, p
      integer             ndata, nn, thin, thinnum  
      parameter           (ndata = 100, d = 5, sde1 = 1.0, sde2 = 1.0)
!cc                        sample size, the number of B-spline basis functions 

      parameter           (thin = 25, nn = 200)
      parameter           (itr = thin*nn, burn = 500)

      double precision    B(ndata), E(ndata), AA(ndata), A, sde
      double precision    DXX(ndata,ndata), XX(ndata,2), zz(ndata)   
      double precision    x_left,  x_cnter, x_right, xma, xint, u1,u2
      double precision    slope_left, inter_left, slope_right, inter_right
      double precision    con1, con2, com_minus, v_left, v_right, bt0, tmk
      double precision    tmean(ndata), tq, lcc(ndata), lev, fix(d), tsum
      double precision    tmp0(ndata), tmp1, tmp2, gamma 
      double precision    gam0(ndata),gam1(ndata), rt, rtmp
      double precision    gam00(ndata), gam11(ndata), gtmp0, gtmp1
      double precision    sgamma, ADXX(ndata, ndata), AXX(ndata, 2)
      double precision    Azz(ndata)

      integer             vv, crt, ip, idx   
      double precision    R, c
      double precision    pp(99), ptmp, mintmp, lk(99), ssp, lksum  
  
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc     I/O part 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cc
!cc Input files 
!cc

	open (unit=10, file='simul_T_M2_N100C1.txt')
	open (unit=11, file='simul_delta_M2_N100C1.txt')
	open (unit=12, file='simul_Z_M2_N100C1.txt')
	open (unit=14, file='simul_bspl_M2_N100C1.txt' )  

!cc     T_D: time, delta: censoring indicator, Z: group variable (0/1), 
!cc     bspl_D: basis expansion of t 

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cc
!cc Output files 
!cc

	open (unit=100, file='beta_bspl_D.txt')
	open (unit=120, file='gamma_D.txt' ) 


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc     Read data 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       call random_seed (put=seed)

      do 600 i = 1, 99 
        pp(i) = i / 100.0 
600   continue     

      do 9 kk = 1, 1
        ip = ndata * (kk-1)

        do 10 j = 1, d	    

          do 11 i = 1, ndata  
 	    read(14, *) ADXX(ip + i, j)       
	    if ( j .eq. 1) then 
              E(i) = 1        
              read (10, *) AXX(ip + i, 1)      
	      read (11, *) AXX(ip + i, 2)  
              read (12, *) Azz(ip + i)  
             endif

11         continue        

10      continue

9     continue
	
      do 2000 llp  = 1, 1      

! ccccccccccccc

        B = 0.0 
	tmean =0.0
	gamma =0.0
        thinnum = 0.0 
	sgamma = 0.0 
        spp = 0.0             
 
       do 50 j = 1, d 

 	 do 51 i = 1, ndata  

	    DXX(i, j) = ADXX(ndata * (llp - 1) + i, j)     
	    if ( j .eq. 1) then 
              XX(i, 1) = AXX(ndata * (llp - 1) + i, 1)  
              XX(i, 2) = AXX(ndata * (llp - 1) + i, 2)
              zz(i) = Azz(ip + i)
            endif

51      continue        

50     continue


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       do 900 ii = 1, itr + burn    
   	
       	 do 500 jj = 1, d
	   c = 0.0 
	   sde = sde2
	   if ( jj .eq. 1 ) sde = sde1  

           do 300 k = 1, ndata
              c = c + ( zz(k) * DXX(k, jj) ) * XX(k, 2)

300        continue	

!c       tmean(jj) = B(jj) / ii + tmean(jj) * (ii - 1) / ii        
!c       tq = tmean(jj)
         tq = 0.0        
         loc = jj 
	  
         if (gamma .eq. 1) then 
  	   call Rgen(tq, loc, B, sde, XX, DXX, zz, E, c, ndata, d, R, num)  
         endif

         if (gamma .eq. 0) then
            if ( jj .eq. 0 ) then  
               call Rgen(tq, 1, B, sde, XX, DXX, zz, E, c, ndata, 1, R, num)  
            else
               R = random_normal()
               R = sde * R  
            endif
	 endif	    
 
         sm = sm + num          
         B(jj) = R      
   
500    continue 

       tsum = 0.0
       thinnum = thinnum + 1

       if ( (thinnum .ge. thin).and.(ii .gt. burn) ) then  
         do 3 kk = 1, d 	
            if ( kk .eq. 1 ) B(kk) = B(kk) + tsum 
            if ( kk .ne. 1 ) B(kk) = B(kk) + tsum         
            if((ii.gt.burn).and.(sh.eq.1)) then
     	      write(100, 9990) ii-burn, kk, B(kk)
            endif   
            if((ii.gt.burn).and.(sh.eq.2)) then
      	      write(200, 9990) ii - burn, kk, B(kk) 
            endif 

3        continue  
       endif 

       if ( gamma .eq. 1 ) call gval_atom(B, DXX, XX, zz, ndata, d, AA)          
       if ( gamma .eq. 0 ) call gval_atom(B, DXX, XX, zz, ndata, 1, AA)  

       do 700 kk = 1, ndata 	 
         tmk = random_exponential()
         E(kk) = tmk / AA(kk)   	    

700    continue 

          gam0 = 0.0
          gam1 = 0.0
          gam00 = 0.0
	  gam11 = 0.0  
	  lk = 0.0 

        do 601 i = 1, 99

           if (gamma.eq.0.0) then
              call gval_atom(B, DXX, XX, zz, ndata, 1, gam00)
              call get_beta(B(1), 1, B, DXX, ndata, 1, gam0)  
           do 602 kk=1, ndata
              lk(i) = lk(i) + ( gam0(kk) * zz(kk) - E(kk) * gam00(kk)) * XX(kk, 2)               

602        continue

           do 603 kk = 1, 1 
              lk(i) = lk(i) - 0.5 * B(1)**2 / (sde**2)  

603        continue
              lk(i) = lk(i) + dlog( 1 - pp(i) ) 

            endif 

           if (gamma.eq.1.0) then
              call gval_atom(B, DXX, XX, zz, ndata, d, gam11)
              call get_beta(B(1), 1, B, DXX, ndata, 1, gam1)  
            do 604 kk=1, ndata
              lk(i) = lk(i) + (gam1(kk) * zz(kk) - E(kk) * gam11(kk)) * XX(kk, 2)               

604         continue

            do 605 kk=1,d 
              lk(i) = lk(i) - 0.5 * B(kk)**2 / (sde**2)  

605         continue

              lk(i) = lk(i) + dlog(pp(i)) 
            endif             

601      continue 

          mintmp = lk(1)
          ptmp = 0.0 

          do 606 i=2,99            
            if(lk(i).lt.mintmp) mintmp = lk(i)          

606       continue

          do 607 i = 1, 99
            lk(i) = lk(i) - mintmp
            lk(i) = exp( lk(i) ) 
            ptmp = ptmp + lk(i)

607       continue

            CALL RANDOM_NUMBER(mintmp) 
            idx = 1 
   	    lksum = lk(1)

          do while (lksum / ptmp.lt.mintmp)
            idx = idx + 1              
	    lksum = lksum + lk(idx) 
          end do 	         

           p = pp(idx)
           if(p.ge.1.0) p = 0.99
	   if(p.le.0.0) p = 0.01		             

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 

          gam0 = 0.0
          gam1 = 0.0
          gam00 = 0.0
	  gam11 = 0.0

	  gtmp0 = 0.0
	  gtmp1 = 0.0

          call gval_atom(B, DXX, XX, zz, ndata, 1, gam00)                  
          call gval_atom(B, DXX, XX, zz, ndata, d, gam11)
	     	  
    	  call get_beta(B(1), 1, B, DXX, ndata, 1, gam0) 
    	  call get_beta(B(1), 1, B, DXX, ndata, d, gam1)
  
        do 800 kk = 1, ndata
	   gtmp0 = gtmp0+(gam0(kk) * zz(kk) - E(kk) * gam00(kk)) * XX(kk,2) 
	   gtmp1 = gtmp1+(gam1(kk) * zz(kk) - E(kk) * gam11(kk)) * XX(kk,2)     	  

800     continue	     

        rtmp = gtmp0 - gtmp1     
        gamma = 0.0

        if ( abs(rtmp) .le. 700 ) then
          rt = p / (dexp(rtmp) * (1 - p) + p) 
          CALL RANDOM_NUMBER(mintmp) 
          if ( mintmp .le. rt ) then
            gamma =1.0 
          endif 	   
        endif 

          if ( rtmp .gt.  700 ) gamma = 0.0 
	  if ( rtmp .le. -700 ) gamma = 1.0  

       if ( (thinnum .ge. thin).and.(ii .gt. burn) ) then  
          if(ii.gt.burn) write(120, 9993) ii - burn, gamma	    
          sgamma = sgamma + gamma / nn
          spp = spp + p / nn                     
          thinnum = 0.0
       endif 

       print *, itr, ii - burn, gamma 
           
       sm = 0.0    
        
900   continue
     
       if ( sh .eq. 1) write(1000, *) llp, sgamma, spp 

2000  continue 

9990  format ( i8, 1x, i4, 1x, f10.5)      
9991  format ( i8, 1x, i4, 1x, f10.5)      
9992  format ( f10.5, 1x, f10.5 )
9993  format ( i8, 1x, f10.5 )      
      end 

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Sub 1
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine Rgen(bt0, loc, B, sde, XX, DXX, zz, E, c, ndata, d, R, num) 

         integer            d, num, loc, ndata
	 double precision   bt0, sde, R, rndata
	 double precision   B(ndata)   
	 double precision   bt, XX(ndata, 2), E(ndata), AA(ndata), A   
	 double precision   DXX(ndata,ndata), zz(ndata)
         double precision   x_left,  x_center1, x_center2, x_right, xma, xint
	 double precision   slope_left, inter_left, slope_right, inter_right
	 double precision   con1, con2, con3, com_minus, v_left, v_right  
         double precision   dd, up, pp, tmp, vv, crt, u1,u2, cbb, c, xmaval
	 
	 parameter(cbb = 0.1)          

         rndata = ndata * 1.0
         v_left = 1 / rndata**0.50 
 	 v_right = 1 / rndata**0.50
      	 
         call getmax_posterior (bt0, loc, B, sde, XX, DXX, zz, E, c, ndata, d, bt) 	   
         xma = bt  
         call getpts(xma, loc,B, v_left, v_right, sde, XX, DXX, zz, E, c, ndata, d, x_left, x_center1, x_center2, x_right, &
                     slope_left, inter_left, slope_right, inter_right, con1, con2, con3, com_minus, xmaval)
         vv = 1 

 	 do 20 while (vv .gt. 0 ) 
	 if ( slope_left .le. cbb )  then  
           v_left = v_left * 1.5 
           call getpts(xma, loc,B, v_left, v_right, sde, XX, DXX, zz, E, c, ndata, d, x_left, x_center1, x_center2, x_right, & 
                       slope_left, inter_left, slope_right, inter_right, con1, con2, con3, com_minus, xmaval)
 	 endif

	 if ( slope_right .ge. -cbb ) then  
           v_right = v_right * 1.5
           call getpts(xma, loc,B, v_left, v_right, sde, XX, DXX, zz, E, c, ndata, d, x_left, x_center1, x_center2, x_right, & 
                       slope_left, inter_left, slope_right, inter_right, con1, con2, con3, com_minus, xmaval)             
	 endif
         if (  abs(slope_left * slope_right) .ge. cbb*cbb ) vv = 0 

20       continue   

         tmp=con1 + con2 + con3   
 	 con1 = con1 / tmp 
	 con2 = con2 / tmp
	 con3 = con3 / tmp 

         crt = 1
 	 num = 0

	 do 30 while (crt .gt. 0 ) 	

         CALL RANDOM_NUMBER(u1)   
         CALL RANDOM_NUMBER(u2)   

        if ( u1 .le. con1 ) then 
             R = ( dlog( u1 * tmp * slope_left ) - inter_left ) / slope_left 
        endif
      
	if ( (u1 .ge. con1).and.(u1 .le. (con1 + con2) ) ) then 
             R = tmp * (u1 - con1) / dexp(xmaval) + x_center1	  
	endif   

	if ( (1 - u1) .le. con3 ) then 
             R = ( dlog( (1 - u1) * tmp * abs(slope_right) ) - inter_right ) / slope_right     
	endif   
        up = xmaval

        if ( R .le. x_center1 ) up = slope_left * R + inter_left 
        if ( R .gt. x_center2 ) up = slope_right * R + inter_right 	 

        A = 0.0	
        call gval(R, loc, B, 0, XX, DXX, zz, E, ndata, d, A)	 
        pp = c * R - A - 0.5 * R**2 / sde**2 - com_minus   

        A = 0.0
        A = dlog(u2)
        num = num + 1
	if ( A .le. (pp - up) ) crt =0       
	
30      continue      

        return
        end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Sub 2
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine getmax_posterior(bt0, loc, B, sde, XX, DXX, zz, E, c, ndata, d, bt)  
         
         integer            loc, d, ndata 
         double precision   bt0, bt, rnadta 
	 double precision   B(ndata) 
	 double precision   XX(ndata,2), DXX(ndata,ndata), zz(ndata), E(ndata)   
	 double precision   A1, A2
	 double precision    crit, lp, llp, tmp, sde, c  
	 
         crit =1
         bt = bt0 
         tmp = c  
         rndata = ndata 
      
         A1 = 0.0
         A2 = 0.0 

	 do 10 while ( abs(crit) .ge. 0.001 / sqrt(rndata) )  
 
         call gval(bt, loc, B, 1, XX, DXX, zz, E, ndata, d, A1)
     	 lp = tmp - A1 - bt / sde**2  

      	 call gval(bt, loc, B, 2, XX, DXX, zz, E, ndata, d, A2)            
         llp = - A2 - 1 / sde**2
        
         crit = lp / llp
	 bt = bt - crit  

10       continue  
         return
         end   

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Sub 3
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine getpts(xma, loc, B, v_left, v_right, sde, XX, DXX, zz, E, c, ndata, d, x_left, x_center1, x_center2, x_right, & 
                  slope_left, inter_left, slope_right, inter_right, con1, con2, con3, com_minus, xmaval) 	 

        integer            d, loc, ndata 
        double precision   B(ndata)
        double precision   xma, v_left, v_right, XX(ndata,2), E(ndata)
        double precision   DXX(ndata,ndata), zz(ndata)
        double precision   sde, x_left, x_center1, x_center2, x_right
        double precision   slope_left, inter_left, slope_right, inter_right
        double precision   con1, con2, con3, com_minus, A0, A1, tmp, xmaval
        double precision   xint, c   

        x_left = xma  - v_left 
        x_right = xma + v_right

        slope_left = 0
        inter_left = 0

	tmp = c 
	A0 = 0 
	A1 = 0 	 

        call gval(x_left, loc, B, 1, XX, DXX, zz, E, ndata, d, A1) 
        call gval(x_left, loc, B, 0, XX, DXX, zz, E, ndata, d, A0) 
	
        slope_left = tmp - A1 - x_left / sde**2
        inter_left = x_left * tmp - A0 - 0.5 * x_left**2 / sde**2
	inter_left = inter_left - slope_left * x_left 

	A0 =0 
	A1 =0 	
        call gval(x_right, loc, B, 0, XX, DXX, zz, E, ndata, d, A0) 
        call gval(x_right, loc, B, 1, XX, DXX, zz, E, ndata, d, A1) 

        slope_right = tmp - A1 - x_right / sde**2
        inter_right = x_right * tmp - A0 - 0.5 * x_right**2 / sde**2
        inter_right = inter_right - slope_right * x_right 

	A0 =0 	
        call gval(xma, loc, B, 0, XX, DXX, zz, E, ndata, d, A0)       		 
        xmaval =  tmp * xma - A0 - 0.5 * xma**2 / sde**2	

	x_center1 = - (inter_left-inter_right) / (slope_left-slope_right)
	x_center2 = - (inter_left-inter_right) / (slope_left-slope_right)
  
        x_center1 =  (xmaval - inter_left) / slope_left
        x_center2 =  (xmaval - inter_right) / slope_right 
 
        com_minus = inter_right 
	if ( inter_left .le. inter_right ) com_minus = inter_left 
	if ( xmaval .le. com_minus ) com_minus = xmaval

        inter_left = inter_left - com_minus 
	inter_right = inter_right - com_minus
	xmaval = xmaval - com_minus

	if ( inter_left .gt. 700 ) inter_left = 700
	if ( inter_left .le. -700 ) inter_left = -700	
	if ( inter_right .gt. 700 ) inter_left = 700	
	if ( inter_right .le. -700 ) inter_left = -700	
	if ( xmaval .gt. 700 ) inter_left = 700
	if ( xmaval .le. -700 ) inter_left = -700

	con1 = dexp( slope_left * x_center1 + inter_left) / slope_left 
	con2 = dexp( xmaval ) * ( x_center2 - x_center1) 
        con3 = dexp( slope_right * x_center2 + inter_right) / abs(slope_right)	

	return
	end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Sub 4
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine gval(bt, loc, B, od, XX, DXX, zz, E, ndata, d, A) 
 
       integer               od, loc, d, qn, ndata      
       double precision      A, Res(ndata), tmp
       double precision      DXX(ndata,ndata), bt, B(ndata) 
       double precision      XX(ndata,2), zz(ndata), E(ndata)
      	 
       A = 0.0              
       ip = 0.0
 
       call get_beta(bt, loc, B, DXX, ndata, d, Res)

       do 10 i=1, ndata 	 
          tmp = 0.0         
          do 20 j = i, ndata         
             tmp = tmp + dexp( zz(j) * Res(i) ) * (zz(j) ** od)

20        continue

          A = A + XX(i,2) * E(i) * ( DXX(i,loc)**od ) * tmp 

10      continue  

        return
        end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Sub 5
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine  gval_atom (B, DXX, XX, zz, ndata, d, AA) 
	 
	 integer          od, d, qn, ndata
         double precision B(ndata), XX(ndata,2), E(ndata)         
	 double precision DXX(ndata,ndata), zz(ndata)
	 double precision AA(ndata), tmp, Res(ndata) 

	 AA = 0.0
	 ip = 0.0
	        
	 call get_beta (B(1), loc, B, DXX, ndata, d, Res)	         	    	

	 do 10 i=1, ndata               	    	
	    tmp = 0.0
            do 20 j = i, ndata   
              if (d .ne. 1) tmp = tmp + dexp (zz(j) * Res(i))
              if (d .eq. 1) tmp = tmp + 1  

20          continue

            AA(i) = tmp 

10       continue

	 return
	 end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Sub 6
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine get_beta(bt,loc,B,DXX,ndata,d,Res)    

        integer               d, loc, ndata      
	double precision      Res(ndata), bt
	double precision      DXX(ndata,ndata), B(ndata) 
 	
	Res = 0.0
!	print *, bt 

       do 100 i =1, ndata       
!	  Res(i) = DXX(i,loc)*B(loc) + DXX(i,ndata)	 	     

      do 200 j =1, d 
	if ( j .ne. loc ) Res(i) = Res(i) +  DXX(i,j)*B(j) 
	if ( j .eq. loc ) Res(i) = Res(i) +  DXX(i,j)*bt
 200  continue
       
100   continue

      if ( d .eq.1) Res = 0.0  
    
      return
      end   


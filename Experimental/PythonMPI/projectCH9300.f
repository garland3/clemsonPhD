        subroutine projectMain(iteration, method, nx,ny)    
      
c        parameter (nx=25,ny=25)
        integer iteration, method, nx,ny, ierror
        real*8 beta,err,eps,w,yl,xl,dx,dy
        real*8 tnew(nx,ny),told(nx,ny),wbeta,betaw
        real*8 tn(nx,ny),tn1(nx,ny),tn12(nx,ny)



        xl=1.0d+00
        yl=0.75d+00
        dx=xl/dfloat(nx-1)
        dy=yl/dfloat(ny-1)
        beta=(dx*dx)/(dy*dy)
 
        w=1.0d+00


c        do 103 lll=1,11
           
            w=w+(0.1d+00 *iteration)
            print*, 'Starting  [method,w] = ' , method, w           

            wbeta=w/(2.0d+00*(beta+1.0d+00))
            betaw=2.0d+00*(1.0d+00+beta)
            
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      Method = 1 Jacobi Point                                         C
C             = 2 Gauss Point                                          C
C             = 3 Gauss Line                                           C
C             = 4 Jacobi Line                                          C
C             = 5 Gauss Line  (ADI)                                    C
C             = 6 Jacobi Line (ADI)                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c           do method=1,6
c The method is passed as a parameter

               eps=1.d-14
               err=0.0d+00

c               print*, "Method Number ", method 

               if (method.gt.4) then
                   do j=2,ny-1
                       do i=2,nx-1
                           tn1(i,j)=0.0d+00
                           tn(i,j)=0.0d+00
                           tn12(i,j)=0.0d+00
                       enddo
                   enddo

                   do j=1,ny
                       tn(nx,j)=1.0d+00
                       tn12(nx,j)=1.0d+00
                       tn1(nx,j)=1.0d+00
                       tn(1,j)=1.0d+00
                       tn12(1,j)=1.0d+00
                       tn1(1,j)=1.0d+00
                   enddo

                   do i=1,nx
                       tn(i,1)=1.0d+00
                       tn12(i,1)=1.0d+00
                       tn1(i,1)=1.0d+00
                       tn(i,ny)=3.0d+00
                       tn12(i,ny)=3.0d+00
                       tn1(i,ny)=3.0d+00
                   enddo

                   goto 111
               endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
C the else part of the if statement  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCc    
               do j=2,ny-1
                   do i=2,nx-1
                       told(i,j)=0.0d+00
                       tnew(i,j)=0.0d+00
                   enddo
               enddo

               do j=1,ny
                   told(nx,j)=1.0d+00
                   tnew(nx,j)=1.0d+00
                   told(1,j)=1.0d+00
                   tnew(1,j)=1.0d+00
               enddo

               do i=1,nx
                   told(i,1)=1.0d+00
                   tnew(i,1)=1.0d+00
                   told(i,ny)=3.0d+00
                   tnew(i,ny)=3.0d+00
               enddo

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
C end the if statement on label 111
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCc   
 111  continue

       if (method.eq.1) call JPoint   (nx,ny,beta,wbeta,err,eps,w,dx,dy,
     c                                   tnew,told,iter)
       if (method.eq.2) call GSPoint  (nx,ny,beta,wbeta,err,eps,w,dx,dy,
     c                                   tnew,told,iter)
       if (method.eq.3) call GSLine   (nx,ny,beta,wbeta,err,eps,w,dx,dy,
     c                                   tnew,told,iter)
       if (method.eq.4) call JLine    (nx,ny,beta,wbeta,err,eps,w,dx,dy,
     c                                   tnew,told,iter)
       if (method.eq.5) call GSLineADI(nx,ny,betaw,beta,err,eps,w,dx,dy,
     c                                   tn,tn1,tn12,iter)
       if (method.eq.6) call JLineADI (nx,ny,betaw,beta,err,eps,w,dx,dy,
     c                                  tn,tn1,tn12,iter)

      call tecplot(method,tnew,told,tn,tn12,tn1,dx,dy,nx,ny,iter)

C end the do loop started on line 33
c       enddo

C end the do loop started on line 17
c 103   continue
        return
       end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      Creating Tecplot Files                                          C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       subroutine tecplot(method,tnew,told,tn,tn12,tn1,dx,dy,nx,ny,iter)
       integer nx,ny,method,iter
       real*8 tnew(nx,ny),told(nx,ny)
       real*8 tn(nx,ny),tn1(nx,ny),tn12(nx,ny),dxx,dyy,dx,dy
       character*16 filename1
       if (iter.ge.1000000) goto 1001
       print*, "Tecplot Files "
        

       if (method.eq.1) filename1='JP.plt'
       if (method.eq.2) filename1='GSP.plt'
       if (method.eq.3) filename1='GSL.plt'
       if (method.eq.4) filename1='JL.plt'
       if (method.eq.5) filename1='GSLADI.plt'
       if (method.eq.6) filename1='JLADI.plt'
       
       open(unit=1,file=filename1,status='unknown',access='append')

       if (method.lt.5) then
           write(1,*) 'VARAIBLES = "X","Y","tnew"'
           write(1,*) 'ZONE F=POINT, I=',nx,', J=',ny
           
            do j=1,ny
            dyy=0.0d+00
            dyy=dyy+dfloat(j-1)*dy
            do i=1,nx
            dxx=0.0d+00
            dxx=dxx+dfloat(i-1)*dx
            write(1,*)dxx,dyy,tnew(i,j)
            enddo
            enddo

           close(unit=1)
       else
       
           write(1,*) 'VARAIBLES = "X","Y","tn1"'
           write(1,*) 'ZONE F=POINT, I=',nx,', J=',ny

            do j=1,ny
            dyy=0.0d+00
            dyy=dyy+dfloat(j-1)*dy
            do i=1,nx
            dxx=0.0d+00
            dxx=dxx+dfloat(i-1)*dx
            write(1,*)dxx,dyy,tn1(i,j)
            enddo
            enddo
            
            close(unit=1)
        
        endif
        
       print*, "Ending Tecplot Files "
 1001  return
       end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      Jacobi Point Method
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

       subroutine JPoint(nx,ny,beta,wbeta,err,eps,w,dx,dy,tnew,told,
     c                                                          iter)
       
       integer nx,ny
       integer iter
       real*8 tnew(nx,ny),told(nx,ny)
       real*8 wbeta,beta,err,eps,w,dx,dy
       
       print*, "JPoint method"
       
        do iter=1,10000000
            if(Modulo(iter,1000000).eq.0) then
                print*, "JPoint iteration number = ", iter
            endif
            
            do j=2,ny-1
                do i=2,nx-1
                    tnew(i,j)=(1.0d+00-w)*told(i,j)+wbeta*
     c (told(i+1,j)+told(i-1,j)+beta*(told(i,j+1)+told(i,j-1)))

                enddo
            enddo

            do j=2,ny-1
                do i=2,nx-1
                    err=err+dabs(told(i,j)-tnew(i,j))
                enddo
            enddo

            err=err/dfloat(nx*ny)

            if (err.lt.eps) then
                print*,"eps and error ", eps,err
                goto 501
            endif

            do mm=2,ny-1
                do mn=2,nx-1
                    told(mn,mm)=tnew(mn,mm)
                enddo
            enddo
        
        enddo
        
 501   open(unit=2,file='iter.txt',status='unknown',access='append')

       write(2,*),'Jpoint    ','iter = ',iter,'omega =  ',w
       
       close(unit=2)

        return
        end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC    Gauss Seidel Point Method
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

       subroutine GSPoint(nx,ny,beta,wbeta,err,eps,w,dx,dy,tnew,told
     c                          ,iter)
       integer nx,ny
       real*8 told(nx,ny),tnew(nx,ny)
       real*8 wbeta,beta,err,eps,w,dx,dy
       
      print*, "Gauss Seidel Point method"

        do iter=1,1000000
            if(Modulo(iter,1000001).eq.0) then
                print*, "Gauss Seidel iteration number = ", iter
            endif
        
        do i=2,nx-1
        do j=2,ny-1
      tnew(i,j)=(1.0d+00-w)*told(i,j)+wbeta*
     c (told(i+1,j)+tnew(i-1,j)+beta*(told(i,j+1)+tnew(i,j-1)))
        enddo
        enddo
      do i=2,nx-1
      do j=2,ny-1
      err=err+dabs(told(i,j)-tnew(i,j))
      enddo
      enddo

      err=err/dfloat(nx*ny)

        if (err.lt.eps) then
        print*,eps,err
        goto 502
        endif

        do mm=2,nx-1
        do mn=2,ny-1
        told(mm,mn)=tnew(mm,mn)
        enddo
        enddo
        enddo

 502   open(unit=2,file='iter.txt',status='unknown',access='append')
       write(2,*),'GSpoint    ','iter = ',iter,'omega =  ',w
       close(unit=2)

        return
        end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC    THE GSL METHOD                                                CCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       subroutine GSLine (nx,ny,beta,wbeta,err,eps,w,dx,dy,tnew,told
     c                          ,iter)
       integer nx,ny,iter
       real*8 told(nx,ny),tnew(nx,ny)
       real*8 wbeta,beta,err,eps,w,dx,dy
        print*, "THE GSL METHOD  "

        do iter=1,1000000
        do j=2,ny-1
        do i=2,nx-1
        tnew(i,j)=(1.0d+00-w)*told(i,j)+wbeta*
     c   (tnew(i+1,j)+tnew(i-1,j)+beta*(told(i,j+1)+tnew(i,j-1)))
        enddo
        enddo
        
         do j=2,ny-1
         do i=2,nx-1
         err=err+dabs(told(i,j)-tnew(i,j))
         enddo
         enddo
         err=err/dfloat(nx*ny)

        if (err.lt.eps) then
        print*,eps,err
        goto 503
        endif


        do mm=2,ny-1
        do mn=2,nx-1
        told(mn,mm)=tnew(mn,mm)
        enddo
        enddo
        enddo

 503   open(unit=2,file='iter.txt',status='unknown',access='append')
       write(2,*),'GSLine    ','iter = ',iter,'omega =  ',w
       close(unit=2)
       
       return
       end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC    THE Jacobi Line METHOD                                        CCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       subroutine JLine (nx,ny,beta,wbeta,err,eps,w,dx,dy,tnew,told
     c                          ,iter)
       integer nx,ny,iter
       real*8 told(nx,ny),tnew(nx,ny)
       real*8 wbeta,beta,err,eps,w,dx,dy
         print*, "Jacobi Line METHOD"
        do iter=1,1000000
             if(Modulo(iter,1000000).eq.4) then
                print*, "Jacobi Line  iteration number = ", iter
            endif
        do j=2,ny-1
        do i=2,nx-1
        tnew(i,j)=(1.0d+00-w)*told(i,j)+wbeta*
     c    (tnew(i-1,j)+tnew(i+1,j)+beta*(told(i,j+1)+told(i,j-1)))
        enddo
        enddo

      do j=2,ny-1
      do i=2,nx-1
      err=err+dabs(told(i,j)-tnew(i,j))
      enddo
      enddo

      err=err/dfloat(nx*ny)


        if (err.lt.eps) then
        print*,eps,err
        goto 504
        endif

        do mm=2,ny-1
        do mn=2,nx-1
        told(mn,mm)=tnew(mn,mm)
        enddo
        enddo
        enddo

 504   open(unit=2,file='iter.txt',status='unknown',access='append')
       write(2,*),'JLine    ','iter = ',iter,'omega =  ',w
       close(unit=2)
       return
       end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC    THE GSL METHOD  (ADI)                                         CCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCC First step
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine GSLineADI(nx,ny,betaw,beta,err,eps,w,dx,dy,tn,tn1,tn12
     c                          ,iter)

       integer nx,ny,iter
       real*8 tn12(nx,ny),a(nx),b(nx),tn(nx,ny)
       real*8 c(nx),d(nx),tn1(nx,ny)
       real*8 betaw,beta,err,eps,w,dx,dy
         print*, "GSL METHOD  (ADI)"
        
       do iter=1,1000000
       do j=2,ny-1
       do i=2,nx-1
       a(i)=-w
       b(i)=betaw
       c(i)=-w
       d(i)=betaw*tn(i,j)+w*(beta*(tn12(i,j-1)+tn(i,j+1))-betaw*tn(i,j))
       end do

        a(1)=0.0d+00
        a(nx-1)=0.0d+00

       d(2)=betaw*tn(2,j)+w*(beta*(tn12(2,j-1)+tn(2,j+1))
     c    -betaw*tn(2,j))+w*tn(1,j)

       d(nx-1)=betaw*tn(nx-1,j)+w*(beta*(tn12(nx-1,j-1)+
     c tn(nx-1,j+1))-betaw*tn(nx-1,j))+w*tn(nx,j)     
     
        c(nx-1)=0.0d+00
        do k=3,nx-1
        b(k)=b(k)-(c(k-1)*a(k-1))/b(k-1)
        d(k)=d(k)-(d(k-1)*a(k-1))/b(k-1)
        end do
        tn12(nx-1,j)=d(nx-1)/b(nx-1)
        do m=nx-2,2,-1
        tn12(m,j)=(d(m)-c(m)*tn12(m+1,j))/b(m)
        end do
        end do
	do n=2,nx-1
        do nk=2,ny-1
        a(nk)=-w*beta
        b(nk)=betaw
        c(nk)=-w*beta
        
      d(nk)=betaw*tn12(n,nk)+w*(tn1(n-1,nk)+
     c tn12(n+1,nk)-betaw*tn12(n,nk))
        end do
        a(1)=0.0d+00
        a(nx-1)=0.0d+00

      d(2)=betaw*tn12(n,2)+w*(tn1(n-1,2)+tn12(n+1,2)
     c -betaw*tn12(n,2))+w*beta*tn12(n,1)
      d(nx-1)=betaw*tn12(n,ny-1)+w*(tn1(n-1,ny-1)+
     c tn12(n+1,ny-1)-betaw*tn12(n,ny-1))+w*beta*
     c tn12(n,ny)
        do ii=3,ny-1
        b(ii)=b(ii)-(c(ii-1)*a(ii-1))/b(ii-1)
        d(ii)=d(ii)-(d(ii-1)*a(ii-1))/b(ii-1)
        end do
	tn1(n,ny-1)=d(ny-1)/b(ny-1)
        do ii=ny-2,2,-1
	tn1(n,ii)=(d(ii)-c(ii)*tn1(n,ii+1))/b(ii)
        end do
        do ni=1,ny
        do ki=1,nx
        tn(ki,ni)=tn1(ki,ni)
        end do
        end do
        do kn=1,ny-1
        do km=1,nx
        err=err+dabs(tn12(km,kn)-tn1(km,kn))
        enddo
        enddo

        err=err/dfloat(nx*ny)

        if (err.lt.eps) then
        print*,eps,err
        goto 505
        endif

        end do
        end do

 505   open(unit=2,file='iter.txt',status='unknown',access='append')
       write(2,*),'GSLineADI     ','iter = ',iter,'omega =  ',w
       close(unit=2)
       return
       end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      Jacobi line Method  (ADI)                                       C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC     First step
CCCCCCCCCCCCCCCCCCCCCCCccccCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       subroutine JLineADI(nx,ny,betaw,beta,err,eps,w,dx,dy,tn,tn1,tn12
     c                          ,iter)
        
        integer nx,ny
        real*8 tn12(nx,ny),a(nx),b(nx),tn(nx,ny)
        real*8 c(nx),d(nx),tn1(nx,ny)
        real*8 betaw,beta,err,eps,w,dx,dy
        
        print*, "Jacobi line Method  (ADI)"

        do iter=1,1000000
        if(Modulo(iter,1000000).eq.6) then
                print*, "Jacobi Line (ADI) iteration number = ", iter
            endif
        do j=2,ny-1
        do i=2,nx-1
        a(i)=-w
        b(i)=betaw
        c(i)=-w
      d(i)=betaw*tn(i,j)+w*(beta*(tn(i,j-1)+tn(i,j+1))
     c -betaw*tn(i,j))
        end do
        a(1)=0.0d+00
        a(nx-1)=0.0d+00
        
      d(2)=betaw*tn(2,j)+w*(beta*(tn(2,j-1)+tn(2,j+1))
     c -betaw*tn(2,j))+w*tn(1,j)

      d(nx-1)=betaw*tn(nx-1,j)+w*(beta*(tn(nx-1,j-1)+
     c tn(nx-1,j+1))-betaw*tn(nx-1,j))+w*tn(nx,j)

        c(nx-1)=0.0d+00

        do k=3,nx-1
        b(k)=b(k)-(c(k-1)*a(k-1))/b(k-1)
        d(k)=d(k)-(d(k-1)*a(k-1))/b(k-1)
        end do
        tn12(nx-1,j)=d(nx-1)/b(nx-1)
        do m=nx-2,2,-1
        tn12(m,j)=(d(m)-c(m)*tn12(m+1,j))/b(m)
        end do
        end do
	do n=2,nx-1
        do nk=2,ny-1
        a(nk)=-w*beta
        b(nk)=betaw
        c(nk)=-w*beta
        
      d(nk)=betaw*tn12(n,nk)+w*((tn12(n-1,nk)+
     c tn12(n+1,nk))-betaw*tn12(n,nk))

        end do
        a(1)=0.0d+00
        a(nx-1)=0.0d+00

      d(2)=betaw*tn12(n,2)+w*((tn12(n-1,2)+tn12(n+1,2))
     c -betaw*tn12(n,2))+w*beta*tn12(n,1)
     
      d(nx-1)=betaw*tn12(n,ny-1)+w*((tn12(n-1,ny-1)+
     c tn12(n+1,ny-1))-betaw*tn12(n,ny-1))+w*beta*tn12(n,ny)

        do ii=3,ny-1
        b(ii)=b(ii)-(c(ii-1)*a(ii-1))/b(ii-1)
        d(ii)=d(ii)-(d(ii-1)*a(ii-1))/b(ii-1)
        end do
	tn1(n,ny-1)=d(ny-1)/b(ny-1)
        do ii=ny-2,2,-1
	tn1(n,ii)=(d(ii)-c(ii)*tn1(n,ii+1))/b(ii)
        end do
        do ni=1,ny
        do ki=1,nx
        tn(ki,ni)=tn1(ki,ni)
        end do
        end do
        do nn=1,ny-1
        do mn=1,nx
       err=err+dabs(tn12(mn,nn)-tn1(mn,nn))
        enddo
        enddo

        err=err/dfloat(nx*ny)

        if (err.lt.eps) then
        print*,eps,err
        goto 506
        endif
        
         enddo
         end do

 506   open(unit=2,file='iter.txt',status='unknown',access='append')
       write(2,*),'JLineADI    ','iter = ',iter,'omega =  ',w
       close(unit=2)
       return
       
       end

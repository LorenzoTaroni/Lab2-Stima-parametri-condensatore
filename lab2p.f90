module mod_3



contains



!*******************************************************************



subroutine chiquad(m,q,x,y,dev,ind,k)



implicit none



    real						   :: t

    integer						   :: n, i

    real, intent(in)			   :: m, q

    integer, intent(out)		   :: ind
   
    real, intent(out)  :: k

	real, dimension(:), intent(in) :: x, y, dev



    n = size(x)

    t = 0.0



    do i = 1,n



      t = t + ((y(i) - m*x(i) - q)**2 / dev(i)**2) 



    end do

    k = t

    if (t >= 3.8157 .and. t <= 21.920) then	! non si rifiuta
                      


                                                ! Test di  chi quadro a due code, si sceglie un livello di significatività pari al 5% e 

						! si verifica che il valore dello stimatore rientri o meno nella zona critica, ossia le 

	  ind = 1				! code della funzione di chi quadro. I due estremi nella condizione della dichiarazione 

						! if valgono per chi quadro a sedici gradi di libertà. (6.908,28.845), alpha = 0.05, nu = 16



  	else	                                ! si rifiuta





      ind = 2





    end if

    

end subroutine chiquad



!*******************************************************************



subroutine hist(a,b,n,r,fn,medg,varg)



implicit none



	real							   :: widbin, p, norm

	integer							   :: numbin, n, k, j

	real, intent(in)				   :: a, b

        real, intent(out) :: medg, varg

    character(len=40), intent(out)	   :: fn

    real, dimension(:), intent(in)	   :: r

    real, dimension(:), allocatable :: c



    !print*, "numero di intervalli:"

    !read*, numbin

    numbin = 50

    allocate(c(numbin))



    widbin = (b - a) / numbin

    c      = 0



    do j = 1,n				! Nei due cicli do si inseriscono le componenti del vettore r nei bin adegutai.



      do k = 1,numbin



        if(a + (k-1)*widbin < r(j) .and. r(j) <= a + k*widbin) then



          c(k) = c(k) + 1

          

		end if



      end do



    end do

    
    medg = sum(r) / n

    
    varg =  sqrt(sum((r-medg)**2) / (n-1))   

    
    norm = sum(c)*widbin

    c = c/norm
    

	print*, "nome file:"

    read*, fn



    open(unit = 1000, file = fn)			! Definizione del file contenente i dati dell'istogramma.



    do k = 1,numbin



      p = widbin * (2*k-1) / 2

      

      write(1000,*) a + p, c(k)



    end do

    

end subroutine hist



!*******************************************************************



subroutine minquad(x,y,sigmay,m,q,varm,varq,cov,rho)



implicit none



  	  real						     :: s00,s01,s11,s10,s20

      integer					   	 :: z, n

	  real, intent(out) 			 :: m,q,varm,varq,cov,rho

      real, dimension(:), intent(in) :: x,y,sigmay



      s00 = 0.0

	  s01 = 0.0

	  s11 = 0.0

	  s10 = 0.0

	  s20 = 0.0



	  n =  size(x)



      do z = 1,n



        s00 = s00 + 1.0/(sigmay(z)**2)

		s01 = s01 + y(z)/(sigmay(z)**2)

		s11 = s11 + x(z)*y(z)/(sigmay(z)**2)

		s10 = s10 + x(z)/(sigmay(z)**2)

		s20 = s20 + (x(z)/sigmay(z))**2

        

	  end do



          m = (s00*s11 - s10*s01) / (s00*s20 - (s10**2))				! Parametri stimati in uscita.

	  varm = s00/(s00*s20 - (s10**2))

	  q = (s20*s01 - s10*s11) / (s00*s20 - (s10**2))

	  varq = s20/(s00*s20 - (s10**2))

	  cov  = -s10/(s00*s20-s10**2)

	  rho  = cov/sqrt(varm*varq)



end subroutine minquad



!*******************************************************************



subroutine gaussinc(n,sigma,z)



implicit none



	integer			      :: i

	real,intent(in)		      :: sigma

	integer,intent(in)	      :: n

	real,dimension(100)	      :: r

	real,dimension(:),intent(out) :: z

	real,dimension(:),allocatable :: w



	allocate(w(n))



	do i = 1,n

	  

	  call random_number(r)

      

	  w(i) = sum(r)/100.0					! La media ha distribuzione gaussiana per il teorema del limite centrale.

          z(i) = sigma * (w(i)-0.5) * sqrt(12.0*100.0)

      

	end do

    

end subroutine gaussinc



!*******************************************************************



subroutine montecarlo(x,med,var)



implicit none



    integer			    :: n

    real, intent(out)	            :: var

    real, intent(inout)	            :: med

    real, dimension(:), intent(in)  :: x

    real, dimension(:), allocatable :: y



	n = size(x)

    allocate(y(n))

    

	med	= sum(x) / n

	y	= (med - x)**2

    var = sum(y) / n / (n-1)



end subroutine montecarlo



!*******************************************************************

end module mod_3



!*******************************************************************



program prog_3



use mod_3



implicit none



	real							:: m, q, varm, varq, cov, rho,&

    								   q_quad, rin, rout,         &

                                       minm, maxm, minq, 	      &

                                       maxq, minrc, maxrc, minv0, &

                                       maxv0, dev, m_s,           &

                                       q_s, varm_s, varq_s, cov_s,&

                                       rho_s, rc_s, v0_s, k, medg, varg

	integer							:: i, d, j, cin, cout, ind,   &

    								   acc, ref, ind_s

    real, parameter					:: pi = 3.14159, del = 0.1,   &

    								   m_par = -6.115E-3,		  &

                                       q_par = 1.293

    character(len=40)				:: fnm, fnq, fnrc, fnv0

    real, dimension(4)				:: a_med, a_var

    real, dimension(14)				:: s, u, v, logv, devlog,  &

    								   t_s, v_s, logv_s, devlog_s 

    real, dimension(:), allocatable :: a_m, a_q, a_varm, a_varq,  &

    								   a_cov, a_rho, a_v0, a_rc



! STUDIO DEI DATI RACCOLTI *****************************************



! Minimi quadrati sui dati raccolti ********************************



	print*, "STUDIO DEI DATI RACCOLTI"

	print*, ""    

	print*, "È stato raccolto un set di 14 misure tempo-tensione."



    open(unit = 701, file = "tempi14clean.txt")

    read(701,*) t_s

    

    open(unit = 702, file = "tensioni14clean.txt")

    read(702,*) v_s

    

	dev 	 = del / sqrt(3.0)

    logv_s   = log(v_s)												! Passaggio al logaritmo per la linearizzazione.

    devlog_s = dev / v_s											! Calcolo della deviazione standard dei logaritmi usando la propagazione della varianza.



    call minquad(t_s,logv_s,devlog_s,m_s,q_s,varm_s,varq_s,cov_s,rho_s)

    print*, ""

    print*, ""

    print*,	"MINIMI QUADRATI SUI DATI RACCOLTI"

    print*, ""

    print*, "Coefficiente angolare:", m_s

    print*, "Varianza coefficiente angolare:", varm_s

    print*, "Termine noto:", q_s

    print*, "Varianza termine noto:", varq_s

    print*, "Covarianza:", cov_s

    print*, "Coefficiente di correlazione:", rho_s



! ******************************************************************



    rc_s = -1.0 / m_s			        ! Le due relazioni legano coefficiente angolare e termine noto rispettivamente a tempo 

    v0_s = exp(q_s)			        ! caratteristico e tensione iniziale.



    print*, ""

    print*, "Tempo caratteristico:", rc_s

    print*, "Tensione iniziale:", v0_s



! Chi quadro sui dati raccolti *************************************



    print*, ""

    print*, ""

    print*, "TEST DI CHI QUADRO SUI DATI RACCOLTI"



    call chiquad(m_s,q_s,v_s,logv_s,devlog_s,ind_s,k)

    

    select case(ind_s)



    case(1)														! Ipotesi non rigettabile.

          

      print*, ""

      print*, "Ipotesi non rigettata."



    case(2)														! Ipotesi rigettabile per alpha = 0.05.



      print*, ""

      print*, "Ipotesti rigettata con significatività del 5%."

        

    end select



! ******************************************************************

    

! ******************************************************************



! SIMULAZIONI DEI SET DI MISURE ************************************



	print*, ""

    print*, ""

	print*, "SIMULAZIONI SET DI MISURE"

    print*, ""

    print*, "Sono stati scelti set da 14 misure ciascuno."

	print*, ""

	print*, "Inserire numero di simulazioni:"

    read*, d

    allocate(a_m(d),a_q(d),a_varm(d),a_varq(d),a_cov(d),a_rho(d),a_v0(d),a_rc(d))



    cin  = 0

    cout = 0

    acc	 = 0

    ref  = 0



	do j = 1,d



      !call random_number(r)



	  !r = 300.0*r

      s = exp(q_s + t_s*m_s)



	  call gaussinc(14,dev,u)										! Generazione delle incertezze gaussiane.



	  v      = s + u												! Assegnazione delle incertezze gaussiane.

      logv   = log(v)												! Passaggio al logaritmo per la linearizzazione.

      devlog = dev / v



      call minquad(t_s,logv,devlog,m,q,varm,varq,cov,rho)

      

      a_m(j)    = m

      a_q(j)    = q

      

      q_quad = ((m-m_s)**2/varm + (q-q_s)**2/varq -			  &

      		   2.0*rho * ((m-m_s)*(q-q_s)) /				  &

               sqrt(varm*varq))/(1-rho**2)



	  open(unit = 101, file = "minimi_quadrati.txt")

      open(unit = 102, file = "q_data_in.txt")

      open(unit = 103, file = "q_data_out.txt")

      

	  write(101,*) m, q, varm, varq, cov, rho, q_quad



      if (q_quad <= 1.0) then



        write(102,*) m, q

        cin = cin + 1



      else



        write(103,*) m, q

        cout = cout + 1



      end if

      

    ! Calcolo del chi quadro per un set simulato *******************

    

	  call chiquad(m,q,t_s,logv,devlog,ind,k)



      select case(ind)



      case(1)														! Ipotesi non rigettabile.

          

        acc = acc + 1



      case(2)														! Ipotesi  rigettabile per alpha = 0.05.



		ref = ref + 1

        

      end select

    

	! **************************************************************



    end do

    

! Generazione degli istogrammi *************************************



	a_rc = -1.0 / a_m												! Le due relazioni legano coefficiente angolare e termine noto rispettivamente a tempo 

    a_v0 = exp(a_q)													! caratteristico e tensione iniziale.



	open(unit = 150, file = "parametri.txt")



	do i = 1,d



      write(150,*) a_m(i), a_q(i), a_rc(i), a_v0(i)



    end do



    minm  = minval(a_m)		! Vengono definiti gli estremi per le ascisse degli istogrammi.

    maxm  = maxval(a_m)

    minq  = minval(a_q)

    maxq  = maxval(a_q)

    minrc = minval(a_rc)

    maxrc = maxval(a_rc)

    minv0 = minval(a_v0)

    maxv0 = maxval(a_v0)



    print*, ""

    print*, ""

    print*, "GENERAZIONE DEGLI ISTOGRAMMI"

    print*, ""

    print*, "Generazione istogramma coefficiente angolare:"

    call hist(minm,maxm,d,a_m,fnm,medg,varg)

    open(unit = 201, file = "file_m")

    write(201,*) "set xrange [", minm, ":", maxm, "]"

    write(201,*) "mu =", medg

    write(201,*) "sigma =", varg

    write(201,*) "gauss(x)=1/(sigma*sqrt(2.*pi))*exp(-(x-mu)**2./(2.*sigma**2))"

    write(201,*) "set arrow from ", medg,",0 to ", medg,",100000 nohead lc rgb 'blue'"

    write(201,*) "plot '", trim(fnm), "' w boxes, gauss(x)"

    print*, ""

    print*, "File da caricare su gnuplot tramite 'load': file_m"



    print*, ""

    print*, "Generazione istogramma termine noto:"

    call hist(minq,maxq,d,a_q,fnq,medg,varg)

    open(unit = 202, file = "file_q")

    write(202,*) "set xrange [", minq, ":", maxq, "]"

    write(202,*) "mu =", medg

    write(202,*) "sigma =", varg

    write(202,*) "gauss(x)=1/(sigma*sqrt(2.*pi))*exp(-(x-mu)**2./(2.*sigma**2))"

    write(202,*) "set arrow from ", medg,",0 to ", medg,",100000 nohead lc rgb 'blue'"

    write(202,*) "plot '", trim(fnq), "' w boxes, gauss(x)"

    print*, ""

    print*, "File da caricare su gnuplot tramite 'load': file_q"



    print*, ""

    print*, "Generazione istogramma tempo caratteristico:"

    call hist(minrc,maxrc,d,a_rc,fnrc,medg,varg)

    open(unit = 203, file = "file_rc")

    write(203,*) "set xrange [", minrc, ":", maxrc, "]"

    write(203,*) "mu =", medg

    write(203,*) "sigma =", varg

    write(203,*) "gauss(x)=1/(sigma*sqrt(2.*pi))*exp(-(x-mu)**2./(2.*sigma**2))"

    write(203,*) "set arrow from ", medg,",0 to ", medg,",100000 nohead lc rgb 'blue'"

    write(203,*) "plot '", trim(fnrc), "' w boxes, gauss(x)"

    print*, ""

    print*, "File da caricare su gnuplot tramite 'load': file_rc"



    print*, ""

    print*, "Generazione istogramma tensione iniziale:"

    call hist(minv0,maxv0,d,a_v0,fnv0,medg,varg)

    open(unit = 204, file = "file_v0")

    write(204,*) "set xrange [", minv0, ":", maxv0, "]"

    write(204,*) "mu =", medg

    write(204,*) "sigma =", varg

    write(204,*) "gauss(x)=1/(sigma*sqrt(2.*pi))*exp(-(x-mu)**2./(2.*sigma**2))"

    write(204,*) "set arrow from ", medg,",0 to ", medg,",100000 nohead lc rgb 'blue'"

    write(204,*) "plot '", trim(fnv0), "' w boxes, gauss(x)"

    print*, ""

    print*, "File da caricare su gnuplot tramite 'load': file_v0"



! ******************************************************************



! Totale dei test di ipotesi sui set simulati **********************



    print*, ""

    print*, ""

    print*, "TEST DI IPOTESI CHI QUADRO SUI DATI SIMULATI"

    print*, ""

    print*, "Set rigettati con significatività del 5%:", ref

    print*, "Percentuale set rigettati:", real(ref)/real(d)*100.0

            

    rin  = real(cin)  / real(d)

    rout = real(cout) / real(d)



    print*, ""

    print*, "Percentuale valori m-q tali che Q^2 <= 1:", rin*100

    print*, ""

    print*, "Percentuale valori m-q tali che Q^2 > 1:", rout*100



! ******************************************************************



! Metodo Montecarlo ************************************************



    call montecarlo(a_m,a_med(1),a_var(1))			! Vengono stimate media e varianza di coefficiente angolare e termine noto di tutte le 

    call montecarlo(a_q,a_med(2),a_var(2))			! simulazioni effettuate, e alla stessa maniera, conoscendo le relazioni che legano il 

    call montecarlo(a_v0,a_med(3),a_var(3))			! coefficiente angolare e il termine noto ripettivamente al tempo caratteristico e alla  

    call montecarlo(a_rc,a_med(4),a_var(4))			! tensione iniziale, anche di quest'ultime si possono calcolare media e varianza delle 

								! stime.

    print*, ""

    print*, ""

    print*, "METODO MONTECARLO"

    print*, ""

    print*, "Media dell stime, coefficiente angolare:", a_med(1)

    print*, "Varianza della media, coefficiente angolare:", a_var(1)

    print*, ""

    print*, "Media delle stime, termine noto:", a_med(2)

    print*, "Varianza della media, termine noto:", a_var(2)

    print*, ""

    print*, "Media delle stime, tensione iniziale:", a_med(3)

    print*, "Varianza della media, tensione iniziale:", a_var(3)

    print*, ""

    print*, "Media delle stime, tempo caratteristico:", a_med(4)

    print*, "Varianza della media, tempo caratteristico:", a_var(4)



! ******************************************************************



end program prog_3

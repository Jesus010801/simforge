      PROGRAM MoveMemb
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      parameter(maxnat=1000000)
      parameter(maxn=100000)
      real*8 num,coorx,coory,coorz
      real*8 coorxn,cooryn,coorzn
      real*8 mayx,menx,mayy,meny,mayz,menz,addx,addy,addz
      character*72 line,car1,car2
      character*20 lin1(maxn)
      integer numat

!______________________________________________________________________
!                  ARCHIVOS USADOS 

      write(6,*)'Cambio en coordenada x (en nm): '
      read(5,*)addx
      write(6,*)'Cambio en coordenada y (en nm): '
      read(5,*)addy      
      write(6,*)'Cambio en coordenada z (en nm): '
      read(5,*)addz      
      open(11,file='dppc512_whole.gro')
      rewind 11
      open(12,file='dppc512_nb.gro')
      rewind 12

      menx=20.0
      mayx=-20.0
      meny=20.0
      mayy=-200.0
      menz=20.0
      mayz=-20.0
      read(11,*)line
      write(12,*)line
      read(11,20)numat 
      write(12,20)numat


!______________________________________________________________________
!        COSTRUCCIÓN DEL NUEVO SISTEMA DE COORDENADAS X,Y,Z 

      do i=1,numat
        read(11,22)car1,coorx,coory,coorz,car2
        coorxn=coorx+addx
        cooryn=coory+addy
        coorzn=coorz+addz
        write(12,22)car1,coorxn,cooryn,coorzn,car2
      enddo
      read(11,24)boxx,boxy,boxz
      write(12,24)boxx,boxy,boxz


20    format(i8)
21    format(a20,3f8.3)
22    format(a20,3f8.3,a25)
23    format(a40)
24    format(3f10.5)



      return
      end

      PROGRAM system
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      parameter(maxnat=1000000)
      parameter(maxn=100000)
      real*8 coorx,coory,coorz,coorxn,cooryn,coorzn
      real*8 boxx,boxy,boxz
      character*72 line3,lin2,line,line2
      character*1 cadn
      integer cad,numat,num1,num2,i,j,k


!______________________________________________________________________

      open(11,file='system_shrink1_em.gro')
      rewind 11
      open(12,file='system_shrink2_em.gro')
      rewind 12
      open(13,file='system_shrink3_em.gro')
      rewind 13
      open(14,file='system_shrink4_em.gro')
      rewind 14
      open(15,file='system_shrink5_em.gro')
      rewind 15
      open(16,file='system_shrink6_em.gro')
      rewind 16
      open(17,file='system_shrink7_em.gro')
      rewind 17
      open(18,file='system_shrink8_em.gro')
      rewind 18
      open(19,file='system_shrink9_em.gro')
      rewind 19
      open(20,file='system_shrink10_em.gro')
      rewind 20
      open(21,file='system_shrink11_em.gro')
      rewind 21
      open(22,file='system_shrink12_em.gro')
      rewind 22
      open(23,file='system_shrink13_em.gro')
      rewind 23
      open(24,file='system_shrink14_em.gro')
      rewind 24
      open(25,file='system_shrink15_em.gro')
      rewind 25
      open(26,file='system_shrink16_em.gro')
      rewind 26
      open(27,file='system_shrink17_em.gro')
      rewind 27
      open(28,file='system_shrink18_em.gro')
      rewind 28
      open(29,file='system_shrink19_em.gro')
      rewind 29
      open(30,file='system_shrink20_em.gro')
      rewind 30
      open(31,file='system_shrink21_em.gro')
      rewind 31
      open(32,file='system_shrink22_em.gro')
      rewind 32
      open(33,file='system_shrink23_em.gro')
      rewind 33
      open(34,file='system_shrink24_em.gro')
      rewind 34
      open(35,file='system_shrink25_em.gro')
      rewind 35
      open(36,file='system_shrink26_em.gro')
      rewind 36
      open(40,file='movie.gro')
      rewind 40

      do i=1,6441
        read(11,18)line
        write(40,18)line
      enddo

      do i=1,6441
        read(12,18)line
        write(40,18)line
      enddo

      do i=1,6441
        read(13,18)line
        write(40,18)line
      enddo

      do i=1,6441
        read(14,18)line
        write(40,18)line
      enddo

      do i=1,6441
        read(15,18)line
        write(40,18)line
      enddo

      do i=1,6441
        read(16,18)line
        write(40,18)line
      enddo

      do i=1,6441
        read(17,18)line
        write(40,18)line
      enddo

      do i=1,6441
        read(18,18)line
        write(40,18)line
      enddo

      do i=1,6441
        read(19,18)line
        write(40,18)line
      enddo

      do i=1,6441
        read(20,18)line
        write(40,18)line
      enddo

      do i=1,6441
        read(21,18)line
        write(40,18)line
      enddo

      do i=1,6441
        read(22,18)line
        write(40,18)line
      enddo

      do i=1,6441
        read(23,18)line
        write(40,18)line
      enddo

      do i=1,6441
        read(24,18)line
        write(40,18)line
      enddo

      do i=1,6441
        read(25,18)line
        write(40,18)line
      enddo

      do i=1,6441
        read(26,18)line
        write(40,18)line
      enddo

      do i=1,6441
        read(27,18)line
        write(40,18)line
      enddo

      do i=1,6441
        read(28,18)line
        write(40,18)line
      enddo

      do i=1,6441
        read(29,18)line
        write(40,18)line
      enddo

      do i=1,6441
        read(30,18)line
        write(40,18)line
      enddo

      do i=1,6441
        read(31,18)line
        write(40,18)line
      enddo

      do i=1,6441
        read(32,18)line
        write(40,18)line
      enddo

      do i=1,6441
        read(33,18)line
        write(40,18)line
      enddo

      do i=1,6441
        read(34,18)line
        write(40,18)line
      enddo

      do i=1,6441
        read(35,18)line
        write(40,18)line
      enddo

      do i=1,6441
        read(36,18)line
        write(40,18)line
      enddo



18    format(a72)
20    format(a15,i5,a24)


      return
      end

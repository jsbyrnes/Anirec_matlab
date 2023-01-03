      subroutine sacread(name,a,ierr)
c  read a SAC-format file
      real*4 ahead
      character*80 name
      character*4 chead(158)
      common/header/ahead(158)
      dimension a(1)
      dimension iah(158)
      equivalence (ahead,iah),(ahead,chead)
      ierr=0
      open(8,file=name,access='direct',recl=512,err=999)
c  read the 158-word header, can use its info 
      read(8,rec=1,err=998)(ahead(i),i=1,128)
c      print *,(ahead(i),i=50,55)
c      print *,(iah(70+i),i=1,10)
      nscan=iah(80)
      irec=2
      read(8,rec=2,err=998)(ahead(i+128),i=1,30),(a(j),j=1,98)
c      print *,'rec=2'
      if(chead(145).ne.' GRN'.and.chead(145).ne.' grn') then
        nword=nscan+158
        nrec=(nword-1)/128+1
        last=nword-128*(nrec-1)
      else       
        print *,'reading greens functions' 
        nword=6*nscan+158
        nrec=(nword-1)/128+1
        last=nword-128*(nrec-1)
      endif
      if(nrec.gt.3) then
        do irec=3,nrec-1
c          print *,'rec=',irec
          read(8,rec=irec,err=998)(a((irec-3)*128+98+i),i=1,128)
        end do
      endif
c      print *,'rec=',nrec
      read(8,rec=nrec,err=998)(a((nrec-3)*128+98+i),i=1,last)
      close(8)
      return
  999 print *,'open error'
      ierr=1
      return
  998 print *,'read error: reading',irec,' out of',nrec
      ierr=1
c      print *,irec
      close(8)
      return
      end
      
      subroutine sacout(name,a)
c  write a SAC-format file
      real*4 ahead
      character*80 name
      character*4 chead(158)
      common/header/ahead(158)
      dimension a(1)
      dimension iah(158)
      equivalence (ahead,iah),(ahead,chead)
      open(8,file=name,access='direct',recl=512)
c  write the 158-word header
      write(8,rec=1)(ahead(i),i=1,128)
c      print *,(ahead(i),i=50,55)
c      print *,(iah(70+i),i=1,10)
      nscan=iah(80)
      if(chead(145).ne.' GRN'.and.chead(145).ne.' grn') then
        nword=nscan+158
        nrec=(nword-1)/128+1
      else       
        print *,'writing greens functions' 
        nword=6*nscan+158
        nrec=(nword-1)/128+1
      endif
      write(8,rec=2)(ahead(i+128),i=1,30),(a(j),j=1,98)
      do irec=3,nrec
        write(8,rec=irec)(a((irec-3)*128+98+i),i=1,128)
      end do
c      print *,nrec,' records written'
      close(8)
      return
      end


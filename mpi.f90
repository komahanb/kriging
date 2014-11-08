!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine MPI_START
        use dimKrig
        implicit none
 !       include 'mpif.h'
        integer :: ierr

!        call MPI_Init(ierr)
!        if (ierr /= MPI_SUCCESS) then
!          stop 'MPI Initializatin error'
!        endif
 !       call MPI_Comm_rank(MPI_COMM_WORLD, id_proc , ierr)
 !       call MPI_Comm_size(MPI_COMM_WORLD, num_proc, ierr)
        
        if(id_proc.eq.0) &
        write(*,'(x,a,i3)')'>> Number of Processors = ',num_proc

!       print *,'Processor ', id_proc+1, 'of',num_proc,'...[OK]'
        end subroutine MPI_START
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine stop_all
        implicit none
!        include 'mpif.h'
        integer :: ierr
 !       call MPI_Barrier(MPI_COMM_WORLD,ierr)
 !       call MPI_Finalize(ierr)
        stop'stop all successfully'
        end subroutine stop_all
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        subroutine spliting(nall,node,nsta,ispl,ispm)
        implicit none
        integer, intent(in) :: nall,node,nsta
        integer,dimension(0:node-1,1000), intent(out) :: ispl
        integer,dimension(0:node-1),      intent(out) :: ispm
        integer :: nbase,nwast,ict,i,j,isum
        ispl  = -1
        ispm  =  0
        nbase = nall / node
        nwast = mod(nall,node)
        if(nbase.ge.998)stop'nbase'
        ict  = 0
        isum = 0
        do i=0,node-1
          do j=1,nbase
             ict = ict + 1
             ispl(i,j) = ict + nsta - 1
             ispm(i)   = ispm(i) + 1
          end do
          if(nwast.ne.0.and.i+1.le.nwast)then
             ict = ict + 1
             ispl(i,nbase+1) = ict + nsta - 1
             ispm(i)         = ispm(i) + 1
          end if
          isum = isum + ispm(i)
!         write(*,'(100i3)')i,(ispl(i,j),j=1,ispm(i))
        end do
        if(ict.ne.nall)stop'spliting'
        if(isum.ne.nall)stop'spliting - ispm'
        end subroutine spliting
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

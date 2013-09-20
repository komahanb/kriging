module timer_mod
  implicit none
  save

  PRIVATE !Everything Private by Default
  
  ! Just the Interface is Public
  PUBLIC :: TimerInit,TimerStart,TimerStop,TimerReport

  ! Maximum # of Timers available
  integer, parameter :: MaxTimer = 10

  ! Starting Clock Time of the Timer
  integer :: InitCount

  ! Current # of Timers in List
  integer :: NCurTimer

  ! Maximum Count of the System Timer
  integer :: MaxSysCount

  ! Timer Names
  character, dimension(MaxTimer) :: TimerName*(80)

  ! Current Time Keepers
  integer, dimension(MaxTimer) :: TStart

  ! Total Time Keepers
  real, dimension(MaxTimer) :: TotalTime

!
!  ------------------------------------------------------------------------------
contains
!  ------------------------------------------------------------------------------
!
  subroutine TimerInit()
    implicit none

    integer :: count,count_rate

    NCurTimer    = 0
    TStart(:)    = 0.0
    TotalTime(:) = 0.0
    call SYSTEM_CLOCK(InitCount,count_rate,MaxSysCount)

  end subroutine TimerInit
!
!  ------------------------------------------------------------------------------
!
  subroutine TimerStart( Name )
    implicit none

    character(*), intent(in) :: Name

    integer :: indx,count_rate,count_max

    ! Find Index
    indx = TimeIndex( Name )

    ! If new timer then add
    if(indx <= 0 )then
      call AddTimer( Name )
      indx = TimeIndex( Name )
    endif

    ! Start Timer last thing
    call SYSTEM_CLOCK(TStart(indx), count_rate, count_max)

    return
  end subroutine TimerStart
!
!  ------------------------------------------------------------------------------
!
  subroutine TimerStop( Name )
    implicit none

    character(*), intent(in) :: Name

    integer :: indx,endcount,count_rate,count_max
    real :: t1

    ! Get Ending Time First Thing
    call SYSTEM_CLOCK(endcount, count_rate, count_max)

    ! Find Index
    indx = TimeIndex( Name )

    if(indx <= 0 )then
     ! call stop_all('TimerStop','Timer does not exist')
    endif

    ! Add to Total
    TotalTime(indx) = TotalTime(indx) +                                    &
                      CalcTime(TStart(indx),endcount,count_rate,count_max)

    return
  end subroutine TimerStop
!
!  ------------------------------------------------------------------------------
!
  subroutine TimerReport()
    implicit none
    
    integer :: i,sl,count,count_rate,count_max
    real :: ttotal
    real, dimension(MaxTimer) :: PercentTime
    integer, dimension(NCurTimer) :: SortIndex

    ! Get the Total Time the Timer ran for
    call SYSTEM_CLOCK(count, count_rate, count_max)

    ttotal = CalcTime(InitCount,count,count_rate,count_max)

    ! Get the Percentage of Total Time
    do i=1,NCurTimer
      PercentTime(i) = 100.0 * (TotalTime(i) / ttotal)
    enddo

    ! Rank them
!!$    if( NCurTimer > 0 )then
!!$      print *,'IndexSort'
!!$      call IndexSort(PercentTime(1:NCurTimer),SortIndex,NCurTimer)
!!$      print *,'Done'
!!$    endif
        
    print *,' '
    print *,'---------------------------------------------------------'
    print *,'                      Timing Report'
    print *,'---------------------------------------------------------'
    print *,'  Timer Name                       Time(sec.)  % of Total'
    print *,'---------------------------------------------------------'
    do i=1,NCurTimer
      !sl = SortIndex( i )
      !write(*,10) TimerName(sl),TotalTime(sl),PercentTime(sl)
       write(*,10) TimerName(i),TotalTime(i),PercentTime(i) 
    enddo
    print *,'---------------------------------------------------------'
    write(*,20) ttotal
    print *,'---------------------------------------------------------'

10  format('   ',A26,' ',F14.3,'       ',F6.2)
20  format('            Total Time = ',F9.3,' seconds')

    return
  end subroutine TimerReport
!
!  ------------------------------------------------------------------------------
!
subroutine AddTimer( Name )
    implicit none

    character(*), intent(in) :: Name

    NCurTimer = NCurTimer + 1

    ! If we have to many timers stop
    if( NCurTimer > MaxTimer )then
      print *,'Increase MaxTimer in Timer Module'
!      call stop_all('AddTimer','MaxTimer is too small')
    endif

    TimerName(NCurTimer) = trim(adjustl(Name))
    TStart(NCurTimer)    = 0
    TotalTime(NCurTimer) = 0.0

    return
  end subroutine AddTimer
!
!  ------------------------------------------------------------------------------
!
  function TimeIndex( Name )
    implicit none

    character(*), intent(in) :: Name
    integer :: TimeIndex
    integer :: i

    ! Set default as not found
    TimeIndex = 0

    do i=1,NCurTimer
      if( trim(adjustl(Name)) == trim(adjustl(TimerName(i))) )then
        TimeIndex = i
        return
      endif
    enddo

    return
  end function TimeIndex
!
!  ------------------------------------------------------------------------------
!
  function CalcTime(StartIndex,EndIndex,CountRate,MaxCount)
    implicit none
    
    integer :: StartIndex,EndIndex,CountRate,MaxCount
    real :: CalcTime
    
    ! If the counter flipped over then we must adjust
    if( EndIndex < StartIndex )then
      CalcTime = real((MaxCount - StartIndex)+EndIndex) / real(CountRate)
    else
      CalcTime = real(EndIndex - StartIndex) / real(CountRate)
    end if
    
    return
  end function CalcTime
!  ------------------------------------------------------------------------------
end module timer_mod
!  ------------------------------------------------------------------------------

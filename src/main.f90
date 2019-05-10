 PROGRAM main

 USE parameters
 USE dts
 USE grid
 USE controls
 USE variables
 USE ibpm
 USE ibpm_pressure
 USE write_responses

 IMPLICIT NONE
 INTEGER :: nfoil, npointtt
 REAL :: prev_lift, prev_drag, variance
 CHARACTER(3) :: file_num
 LOGICAL :: threshold_reached
 !------------------------------- start of program -----------------------!
 ! read parameters
 CALL parameters_input

 ! print stuffs
 WRITE(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
 WRITE(*,*) '   Fast Immersed Boundary Projection Metthod: Forward Run Version        '
 WRITE(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
 WRITE(*,*) '   Developed by Taira & Colonius in 2007                                 '
 WRITE(*,*) '   Updated by Hsieh-Chen Tsai in 2015                                    '
 WRITE(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
 WRITE(*,*) ''

 nfoil = start_n -1
 OPEN(unit=250,file="output/responds/gen_force.dat",form="formatted",status="replace")
 CLOSE(250)

 DO WHILE (nfoil .lt. end_n)
    nfoil = nfoil+1
    WRITE(*,*) '===> Running Forward Simulation'

    ! ********* pre-processing before running a simulation (begins) **********
    ! setting up grid
    WRITE(file_num,"(I3.3)") nfoil
    WRITE(*,*) "===>reading ib"//file_num//".inp"

    OPEN(unit=42,file="input/ib"//file_num//".inp",form='formatted',status='old')
    READ(42,*) npointtt
    CLOSE(42)
    IF (npointtt==0) THEN
        WRITE(*,*) "===> INVALID AIRFOIL-WRITING 0"
        OPEN(unit=250,file="output/responds/gen_force.dat",form="formatted",status="unknown",position="append")
        WRITE(250,*) 0, 0
        CLOSE(250)
        CYCLE
    END IF
    
    CALL grid_setup_eulerian_grid(n_foil = nfoil)

    ! setting up controls
    CALL controls_setup_control

    ! setup variables
    WRITE(*,*) 'Setup variables ....'   ! for debugging
    CALL variables_setup_variables

    ! ********* pre-processing before running a simulation (ends) **********

    ! *************** farward simulation (begins) ******************

    innerloop: DO WHILE ( var%it .le. istop )
    ! Increment time forward
        IF (cfl_flag) THEN
            WRITE(*,*) "===> Force = NaN somehow, writing 0"
            OPEN(unit=250,file="output/responds/gen_force.dat",form="formatted",status="unknown",position="append")
            WRITE(250,*) 0, 0
            CLOSE(250)
            WRITE(*,*) "---------------------STARTING NEW SIMULATION-------------------------"
            cfl_flag = .FALSE.
            exit innerloop
        END IF
        var%it  = var%it  + 1

        it_advance_show = 100
        IF ( MOD(var%it,it_advance_show).eq.0 ) THEN
          WRITE(*,*) "...Advancing to itime =", var%it
        END IF

        !WRITE(*,*) 'Calling ibpm_solver at icpu =', icpu   ! for debugging
        CALL ibpm_solver( dtype_data = var )

        ! responses are output when in fwd and bwd sim, but not opt
        !WRITE(*,*) 'Write output responses at icpu =', icpu ! for debugging
        CALL write_responses_write_file( dtype_data = var, last = .FALSE. )

        IF ((mod(var%it,isave).eq.0).or.(var%it.eq.istop)) THEN  ! save every isave time steps
            CALL variables_write_variables( it = var%it, dtype_data = var )
            ! Compute the pressure if required
            IF (compute_pressure) THEN
                CALL ibpm_pressure_calculate_pressure( dtype_data = var )
                CALL ibpm_pressure_write_pressure( dtype_data = var )
            END IF
        END IF

        ! --------------------Debugging-------------------
        ! WRITE(*,*) "....NOW Drag = ",now_drag
        ! WRITE(*,*) "....NOW Lift = ",now_lift
        ! WRITE(*,*) "....PRE Drag = ",prev_drag
        ! WRITE(*,*) "....PRE Lift = ",prev_lift
        ! WRITE(*,*) "....Drag Diff = ", ABS((prev_drag/now_drag)-1)
        ! WRITE(*,*) "....Lift Diff = ", ABS((prev_lift/now_lift)-1)
        ! WRITE(*,*) ABS((prev_lift/now_lift)-1).lt.threshold*DT
        ! --------------------Debugging-------------------

        IF (use_threshold) THEN  ! save & stop when force change amount under threshold
            threshold_reached = ((ABS((prev_drag/now_drag)-1).lt.threshold*DT).and.(ABS((prev_lift/now_lift)-1).lt.threshold*DT))
            IF (threshold_reached.or.(var%it.eq.istop)) THEN
                
                CALL write_responses_write_file( dtype_data = var, last = .TRUE. )

                CALL variables_write_variables( it = var%it, dtype_data = var )

                ! Compute the pressure if required
                IF (compute_pressure) THEN
                    CALL ibpm_pressure_calculate_pressure( dtype_data = var )
                    CALL ibpm_pressure_write_pressure( dtype_data = var )
                END IF
                WRITE(*,*) "........................END..............."
                EXIT
            END IF
            prev_drag = now_drag
            prev_lift = now_lift
        END IF
    END DO innerloop


    CALL variables_destroy_variables
    CALL controls_destroy_control
    CALL grid_destroy_grid
 END DO
 print *, char(7)
    ! *************** farward simulation (ends) ******************

    WRITE(*,*) 'Finish the run. Destorying variables'
    ! destory grid, control, constraints, and variables
    ! CALL variables_destroy_variables
    ! CALL controls_destroy_control
    ! CALL grid_destroy_grid

    WRITE(*,*) 'Variables destoried'

!------------------------------- end of program -----------------------!

 END PROGRAM main

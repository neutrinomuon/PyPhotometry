! ###########################################################################
! ###  FADO - Fitting Analysis using Differential-evolution Optimisation  ###
! ###########################################################################

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     DIRECTORY : F_FADO -> Main Photometric Filters                        !
!                                                                           !
!     RESUME : This directory contains the subroutines used for             !
!              computing filters. Many of these routines are also used      !
!              in the FADO code.                                            !
!                                                                           !
!              The Transmission curves of these filters may be found        !
!              at the following website:                                    !
!                                                                           !
!                         http://voservices.net/filter/			    !
! 									    !
!              It will be used to compute Photometric magnitudes and	    !
!              colors.							    !
! 									    !
!              We have an observatory dedicated queue which contains	    !
!              the IFU and SLT stare mode OB for each of the standard	    !
!              star listed below. The used slit width is the 5 arcsec	    !
!              one for the 3 arms. The selected readout mode is		    !
!              100kH,1x1,high in both the UVB and VIS arms. The user	    !
!              who is interested in specphot observations with		    !
!              different slit width and readout modes, or with a	    !
!              different strategy than those adopted in the instrument	    !
!              calibration plan (see User Manual), should request the	    !
!              necessary telescope time during phase 1 and provide the	    !
!              relevant OB during phase 2.  The current fully		    !
!              flux-calibrated available spectrophotometric standard	    !
!              stars are marked with a blue background color. Note	    !
!              that BD+17 4708 should not be observed anymore due to	    !
!              its binarity.                                                !
!                                                                           !
!              http://www.eso.org/sci/facilities/paranal/instruments/xshooter/tools/specphot_list.html
!                                                                           !
!              BD+17 4708 == ftp://ftp.stsci.edu/cdbs/current_calspec/bd_17d4708_stisnic_003.ascii 
!                                                                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Jean@Porto - Mon Apr 16 08:21:34 WEST 2009 ++++++++++++++++++++++++++++++++

C     -*- fortran -*-
C     This file is autogenerated with f2py (version:1.24.4)
C     It contains Fortran 77 wrappers to fortran functions.

      subroutine f2pywrapintegralall (integralallf2pywrap, sxvalue
     &s, syvalues, lambda_i, lambda_f, n_lambda, iskeepon, int_type, ver
     &bosity)
      external integralall
      real(kind=8) lambda_i
      real(kind=8) lambda_f
      integer(kind=4) n_lambda
      integer(kind=4) iskeepon
      integer(kind=4) int_type
      integer(kind=4) verbosity
      real(kind=8) sxvalues(n_lambda)
      real(kind=8) syvalues(n_lambda)
      real(kind=8) integralallf2pywrap
      real(kind=8)  integralall
      integralallf2pywrap = integralall(sxvalues, syvalues, lambda
     &_i, lambda_f, n_lambda, iskeepon, int_type, verbosity)
      end


      subroutine f2pywrapslin_interp (slin_interpf2pywrap, xorg_in
     &i, xorg_fin, yorg_ini, yorg_fin, newvalue)
      external slin_interp
      real(kind=8) xorg_ini
      real(kind=8) xorg_fin
      real(kind=8) yorg_ini
      real(kind=8) yorg_fin
      real(kind=8) newvalue
      real(kind=8) slin_interpf2pywrap
      real(kind=8)  slin_interp
      slin_interpf2pywrap = slin_interp(xorg_ini, xorg_fin, yorg_i
     &ni, yorg_fin, newvalue)
      end


      subroutine f2pywrapslog_interp (slog_interpf2pywrap, xorg_in
     &i, xorg_fin, yorg_ini, yorg_fin, newvalue)
      external slog_interp
      real(kind=8) xorg_ini
      real(kind=8) xorg_fin
      real(kind=8) yorg_ini
      real(kind=8) yorg_fin
      real(kind=8) newvalue
      real(kind=8) slog_interpf2pywrap
      real(kind=8)  slog_interp
      slog_interpf2pywrap = slog_interp(xorg_ini, xorg_fin, yorg_i
     &ni, yorg_fin, newvalue)
      end


      subroutine f2pywraplamb_effective (lamb_effectivef2pywrap, t
     &_lambda, t_fluxes, ntlambda, l_lambda, s_fluxes, n_lambda, iskeepo
     &n, int_type)
      external lamb_effective
      integer(kind=4) ntlambda
      integer(kind=4) n_lambda
      integer(kind=4) iskeepon
      integer(kind=4) int_type
      real(kind=8) t_lambda(ntlambda)
      real(kind=8) t_fluxes(ntlambda)
      real(kind=8) l_lambda(n_lambda)
      real(kind=8) s_fluxes(n_lambda)
      real(kind=8) lamb_effectivef2pywrap
      real(kind=8)  lamb_effective
      lamb_effectivef2pywrap = lamb_effective(t_lambda, t_fluxes, 
     &ntlambda, l_lambda, s_fluxes, n_lambda, iskeepon, int_type)
      end


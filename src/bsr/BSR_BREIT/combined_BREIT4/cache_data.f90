!======================================================================
      Subroutine cache_data(read)
!======================================================================
      Use bsr_breit
      Use coef_list,     only: ntrm
      Use term_exp,      only: kt1, kt2, IP_kt1, IP_kt2

      Integer :: read !write/read (0/1) switch

      if(read.eq.0) then

        ntrm_c = ntrm
        kt1_c = kt1
        kt2_c = kt2

        if(allocated(joper_c)) Deallocate(joper_c)
        Allocate(joper_c(noper))
        joper_c = joper
        if(allocated(JT_oper_c)) Deallocate(JT_oper_c)
        Allocate(JT_oper_c(ntrm_c,noper))
        JT_oper_c = JT_oper
        if(allocated(IP_kt1_c)) Deallocate(IP_kt1_c)
        Allocate(IP_kt1_c(kt1_c))
        IP_kt1_c = IP_kt1
        if(allocated(IP_kt2_c)) Deallocate(IP_kt2_c)
        Allocate(IP_kt2_c(kt2_c))
        IP_kt2_c = IP_kt2

      elseif(read.eq.1) then

        ntrm = ntrm_c
        kt1 = kt1_c
        kt2 = kt2_c

        joper = joper_c
        if(allocated(JT_oper)) Deallocate(JT_oper)
        Allocate(JT_oper(ntrm,noper))
        JT_oper = JT_oper_c
        if(allocated(IP_kt1)) Deallocate(IP_kt1)
        Allocate(IP_kt1(kt1))
        IP_kt1 = IP_kt1_c
        if(allocated(IP_kt2)) Deallocate(IP_kt2)
        Allocate(IP_kt2(kt2))
        IP_kt2 = IP_kt2_c

      endif

      End Subroutine cache_data

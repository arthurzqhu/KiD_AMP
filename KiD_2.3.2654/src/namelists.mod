  /  i   k820309    s          18.0        dS^                                                                                                          
       namelists.f90 NAMELISTS                                                           #FILEIN    #FILEOUT    FILEIN FILEOUT                                                                     #IIWARM    #KID_OUTDIR    #KID_OUTFILE    #OVC_FACTOR    #L_REUSE_THOMPSON_LOOKUP 	   IIWARM KID_OUTDIR KID_OUTFILE OVC_FACTOR L_REUSE_THOMPSON_LOOKUP                                                                        
            #L_MPHYS    #L_ADVECT    #L_DIVERGE    #L_PUPDATE    #L_FIX_QV    #L_NOMPHYS_QV    #L_NOADV_QV    #L_POSADV_QV    #L_FIX_THETA    #L_NOMPHYS_THETA    #L_NOADV_THETA    #L_NOADV_HYDROMETEORS    #L_NODIV_HYDROMETEORS    #L_SEDIMENT    #ISURFACE    #L_NOADV_AEROSOLS    #L_NODIV_AEROSOLS    #L_FIX_AEROSOLS    #L_SED_ULT    #L_DIVERGE_ADVECTION    #L_PERIODIC_BOUND    #L_FORCE_POSITIVE     L_MPHYS L_ADVECT L_DIVERGE L_PUPDATE L_FIX_QV L_NOMPHYS_QV L_NOADV_QV L_POSADV_QV L_FIX_THETA L_NOMPHYS_THETA L_NOADV_THETA L_NOADV_HYDROMETEORS L_NODIV_HYDROMETEORS L_SEDIMENT ISURFACE L_NOADV_AEROSOLS L_NODIV_AEROSOLS L_FIX_AEROSOLS L_SED_ULT L_DIVERGE_ADVECTION L_PERIODIC_BOUND L_FORCE_POSITIVE                                                                                                                                                             !            #DT "   #DG_DT #   #MPHYS_SCHEME $   #MPHYS_VAR %   #WCTRL &   #ZCTRL '   #TCTRL (   #PCTRL_Z )   #PCTRL_V *   #PCTRL_T +   #IPCTRL ,   #XCTRL -   #LHF_CTRL .   #SHF_CTRL /   #DIAGLEVEL 0   #DGSTART 1   DT DG_DT MPHYS_SCHEME MPHYS_VAR WCTRL ZCTRL TCTRL PCTRL_Z PCTRL_V PCTRL_T IPCTRL XCTRL LHF_CTRL SHF_CTRL DIAGLEVEL DGSTART                                                                                                                               2            #NUM_H_MOMENTS 3   #NUM_H_BINS 4   #H_SHAPE 5   #MOM_INIT 6   #H_NAMES 7   #MOM_NAMES 8   #MOM_UNITS 9   #NUM_AERO_MOMENTS :   #NUM_AERO_BINS ;   #AERO_MOM_INIT <   #AERO_N_INIT =   #AERO_SIG_INIT >   #AERO_RD_INIT ?   #AERO_NAMES @   #IMOMC1 A   #IMOMC2 B   #IMOMR1 C   #IMOMR2 D   #DONUCLEATION E   #DOCONDENSATION F   #DOCOLLISIONS G   #DOSEDIMENTATION H   #CLOUD_INIT I   #RAIN_INIT J   NUM_H_MOMENTS NUM_H_BINS H_SHAPE MOM_INIT H_NAMES MOM_NAMES MOM_UNITS NUM_AERO_MOMENTS NUM_AERO_BINS AERO_MOM_INIT AERO_N_INIT AERO_SIG_INIT AERO_RD_INIT AERO_NAMES IMOMC1 IMOMC2 IMOMR1 IMOMR2 DONUCLEATION DOCONDENSATION DOCOLLISIONS DOSEDIMENTATION CLOUD_INIT RAIN_INIT                                                                                                                                                                       K            #INPUT_FILE L   #L_INPUT_FILE M   #IFILETYPE N   #ICASE O   INPUT_FILE L_INPUT_FILE IFILETYPE ICASE                                                                    P     
                                                      Q     
                                                      R     
                                                      S     
       MPHYS_ID                                                T     
  
     IMOMC1 IMOMC2 IMOMR1 IMOMR2 DONUCLEATION DOCONDENSATION DOCOLLISIONS DOSEDIMENTATION CLOUD_INIT RAIN_INIT           @                                U     d                 @@                                 A                      @@                                 B                      @@                                 C                      @@                                 D                      @@                                 E                      @@                                 F                      @@                                 G                      @@                                 H                      @@                                 I                   	      p          p            p                                    @@                                 J                   	      p          p            p                                    D@                                 3                         p          p            p                                    D@                                 4                         p          p            p                                    D@                                 5                   	      p          p            p                                    D@                                6                   
      p          p            p                          +          D@                                7            
             p          p            p                                  +          D@                                8            
             p          p            p                                  +          D@                                9            
             p          p            p                                            D@                                 :                         p          p            p                                    D@                                 ;                         p          p            p                                    D@                                <                   
      p          p            p                                    D@                                =                   
      p          p            p                                    D@                                >                   
      p          p            p                                    D@                                ?                   
      p          p            p                          +          D@                                @            
             p          p            p                                            D@                                 "     	                 D@                                #     
                 D@                                $                      D@                                 %                      D@                                &                   
      p          p            p                                    D@                                '                   
      p          p            p                                    D@                                (                   
      p          p            p                                    D@                                )                   
      p          p            p                                    D@                                *                   
      p          p            p                                    D@                                +                   
      p          p            p                                    D@                                 ,                      D@                                -                   
      p          p            p                                    D@                                .                   
      p          p            p                                    D@                                /                   
      p          p            p                                    D@                                 0                      D@                                1     
                 D@ @                              L     d                 D@                                 M                      D@                                 N                      D@                                 O                      D@                                                       D@                                                       D@                                                       D@                                                       D@                                                       D@                                                       D@                                                       D@                                                       D@                                                       D@                                                       D@                                                       D@                                                       D@                                                       D@                                                       D@                                                       D@                                                       D@                                                       D@                                                       D@                                                       D@                                                       D@                                                       D@                                                        D@                                 	                      @                                 V                                                         W                                                      1                                             X                                                      2                                             Y                                                      8                                             Z                                                      5                                             [                                                      7                                             \                                                      6                                             ]                                       	               9                                             ^                                       
               10                                             _                                                      11                                             `                                                      12                                             a                                                      13                                             b                                       �              998          D@                                                       D@                                     �                 D@                                     �                 @@                                    
                 @ @                              c     �                 @ @                              d     �                 @ @                              e     �                 @ @                              f     �                 D@                                     �                 D@ @                                   �                  @                                 g            #         @                                   h                        �          fn#fn    �   r      NAMELISTTOUSE    2  �      ADDCONTROL    &  r     SWITCH    �  �     CONTROL    n  N     MPHYS    �
  �      CASE    t  @   J   PARAMETERS    �  @   J   SWITCHES    �  @   J   SWITCHES_BIN    4  I   J  HEADER_DATA    }  �   J  MICRO_PRM %   '  @       MPHYS_ID+HEADER_DATA !   g  @       IMOMC1+MICRO_PRM !   �  @       IMOMC2+MICRO_PRM !   �  @       IMOMR1+MICRO_PRM !   '  @       IMOMR2+MICRO_PRM '   g  @       DONUCLEATION+MICRO_PRM )   �  @       DOCONDENSATION+MICRO_PRM '   �  @       DOCOLLISIONS+MICRO_PRM *   '  @       DOSEDIMENTATION+MICRO_PRM %   g  �       CLOUD_INIT+MICRO_PRM $   �  �       RAIN_INIT+MICRO_PRM )   �  �       NUM_H_MOMENTS+PARAMETERS &   #  �       NUM_H_BINS+PARAMETERS #   �  �       H_SHAPE+PARAMETERS $   K  �       MOM_INIT+PARAMETERS #   �  �       H_NAMES+PARAMETERS %   {  �       MOM_NAMES+PARAMETERS %     �       MOM_UNITS+PARAMETERS ,   �  �       NUM_AERO_MOMENTS+PARAMETERS )   G  �       NUM_AERO_BINS+PARAMETERS )   �  �       AERO_MOM_INIT+PARAMETERS '   o  �       AERO_N_INIT+PARAMETERS )     �       AERO_SIG_INIT+PARAMETERS (   �  �       AERO_RD_INIT+PARAMETERS &   +  �       AERO_NAMES+PARAMETERS    �  @       DT+PARAMETERS !     @       DG_DT+PARAMETERS &   G  @       MPHYS_SCHEME+SWITCHES #   �  @       MPHYS_VAR+SWITCHES    �  �       WCTRL+SWITCHES    [  �       ZCTRL+SWITCHES    �  �       TCTRL+SWITCHES !   �  �       PCTRL_Z+SWITCHES !     �       PCTRL_V+SWITCHES !   �  �       PCTRL_T+SWITCHES     ?  @       IPCTRL+SWITCHES      �       XCTRL+SWITCHES "     �       LHF_CTRL+SWITCHES "   �  �       SHF_CTRL+SWITCHES %   ;  @       DIAGLEVEL+PARAMETERS #   {  @       DGSTART+PARAMETERS $   �  @       INPUT_FILE+SWITCHES &   �  @       L_INPUT_FILE+SWITCHES #   ;   @       IFILETYPE+SWITCHES    {   @       ICASE+SWITCHES !   �   @       L_MPHYS+SWITCHES "   �   @       L_ADVECT+SWITCHES #   ;!  @       L_DIVERGE+SWITCHES #   {!  @       L_PUPDATE+SWITCHES "   �!  @       L_FIX_QV+SWITCHES &   �!  @       L_NOMPHYS_QV+SWITCHES $   ;"  @       L_NOADV_QV+SWITCHES %   {"  @       L_POSADV_QV+SWITCHES %   �"  @       L_FIX_THETA+SWITCHES )   �"  @       L_NOMPHYS_THETA+SWITCHES '   ;#  @       L_NOADV_THETA+SWITCHES .   {#  @       L_NOADV_HYDROMETEORS+SWITCHES .   �#  @       L_NODIV_HYDROMETEORS+SWITCHES $   �#  @       L_SEDIMENT+SWITCHES "   ;$  @       ISURFACE+SWITCHES *   {$  @       L_NOADV_AEROSOLS+SWITCHES *   �$  @       L_NODIV_AEROSOLS+SWITCHES (   �$  @       L_FIX_AEROSOLS+SWITCHES '   ;%  @       L_SED_ULT+SWITCHES_BIN -   {%  @       L_DIVERGE_ADVECTION+SWITCHES *   �%  @       L_PERIODIC_BOUND+SWITCHES *   �%  @       L_FORCE_POSITIVE+SWITCHES 1   ;&  @       L_REUSE_THOMPSON_LOOKUP+SWITCHES     {&  @       IMPHYS+SWITCHES '   �&  q       IMPHYS_LEM2_4+SWITCHES (   ,'  q       IMPHYS_TAU_BIN+SWITCHES +   �'  q       IMPHYS_THOMPSON09+SWITCHES +   (  q       IMPHYS_THOMPSON06+SWITCHES +   (  q       IMPHYS_THOMPSON07+SWITCHES 0   �(  q       IMPHYS_MORR_TWO_MOMENT+SWITCHES &   a)  q       IMPHYS_UM7_3+SWITCHES %   �)  r       IMPHYS_WSM6+SWITCHES %   D*  r       IMPHYS_WDM6+SWITCHES #   �*  r       IMPHYS_4A+SWITCHES $   (+  r       IMPHYS_AMP+SWITCHES $   �+  s       ITEST_CASE+SWITCHES    ,  @       IIWARM    M,  @       KID_OUTDIR    �,  @       KID_OUTFILE    �,  @       OVC_FACTOR    -  @       FILENAME    M-  @       FILENAMEIN    �-  @       FILENAMEOUT    �-  @       NAMELISTIN    .  @       FILEIN    M.  @       FILEOUT    �.  @       FEXIST    �.  H       READ_NAMELIST 
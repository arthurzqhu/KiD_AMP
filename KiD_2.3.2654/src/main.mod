  S  <   k820309    s          18.0        �|^                                                                                                          
       main.f90 MAIN                                                     
                                                           
                                                           
                @       �                                  
       DT DG_DT NX NZ                                                     
       READ_NAMELIST                                                     
       TIME TIME_STEP N_TIMES                                                     
       READ_PROFILES                                                     
       INTERPOLATE_INPUT INTERPOLATE_FORCING                                                	     
       SAVE_DIAGNOSTICS_1D SAVE_DIAGNOSTICS_2D WRITE_DIAGNOSTICS QUERY_DGSTEP                                                
     
       CALC_DERIVED_FIELDS                                                     
       ADVECT_COLUMN                                                     
       MPHYS_COLUMN                                                     
       STEP_COLUMN                                                     
       DIVERGE_COLUMN                �                                      u #READ_PROFILES_FILE    #READ_PROFILES_STANDARD                      @                                'p                    #H_ID    #NMOMENTS    #NBINS    #MOMENTS                 �                                                               �                                                              �                                                             �                                                          
            &                   &                                                                                            	                                                       
                                                                                                          1                                                                                    x               120#         @                                                                 @                                     
                 @                                                                                                     #         @                                                      #INPUT_TYPE                                                           #         @                                                        #         @                                  !                     #         @                                  "                     #         @                                  #                     #         @                                  $                     #         @                                  %                     #         @                                  &                    #SCHEME_ID '             
                                 '           #         @                                  (                    #SCHEME_ID )             
                                  )           #         @                                  *                     #         @                                  +                                                                 ,                                                        -                      @ @                              .     d                   @                               /                      @ @                               0                                                        1                                                        2                                                        3                        @                               4                             �                            5     �      p              p i        p x         & p         p          p            p x         p          p                          #SPECIES                                                6            #         @                                                     #FILE 7   #INPUT_TYPE 8                                            7                     1                                            8            #         @                                                     #ICASE 9             
                                  9           #         @                                   :                        �         fn#fn    �   @   J   TYPEKIND    �   @   J   SWITCHES !   6  @   J   COLUMN_VARIABLES    v  O   J  PARAMETERS    �  N   J  NAMELISTS      W   J  RUNTIME    j  N   J  SET_PROFILES    �  f   J  INTERPOLATION      �   J  DIAGNOSTICS    �  T   J  DERIVED_FIELDS $   �  N   J  ADVECTION_INTERFACE     G  M   J  MPHYS_INTERFACE    �  L   J  STEPFIELDS    �  O   J  DIVERGENCE /   /  t       gen@READ_PROFILES+SET_PROFILES &   �  �       SPECIES+CLASS_SPECIES +   #  H   a   SPECIES%H_ID+CLASS_SPECIES /   k  H   a   SPECIES%NMOMENTS+CLASS_SPECIES ,   �  H   a   SPECIES%NBINS+CLASS_SPECIES .   �  �   a   SPECIES%MOMENTS+CLASS_SPECIES    �  @       DT+PARAMETERS !   �  @       DG_DT+PARAMETERS    '  q       NX+PARAMETERS    �  s       NZ+PARAMETERS (   	  H       READ_NAMELIST+NAMELISTS    S	  @       TIME+RUNTIME "   �	  @       TIME_STEP+RUNTIME     �	  @       N_TIMES+RUNTIME 0   
  X       INTERPOLATE_INPUT+INTERPOLATION ;   k
  @   a   INTERPOLATE_INPUT%INPUT_TYPE+INTERPOLATION 2   �
  H       INTERPOLATE_FORCING+INTERPOLATION 0   �
  H       SAVE_DIAGNOSTICS_1D+DIAGNOSTICS 0   ;  H       SAVE_DIAGNOSTICS_2D+DIAGNOSTICS .   �  H       WRITE_DIAGNOSTICS+DIAGNOSTICS )   �  H       QUERY_DGSTEP+DIAGNOSTICS 3     H       CALC_DERIVED_FIELDS+DERIVED_FIELDS 2   [  W       ADVECT_COLUMN+ADVECTION_INTERFACE <   �  @   a   ADVECT_COLUMN%SCHEME_ID+ADVECTION_INTERFACE -   �  W       MPHYS_COLUMN+MPHYS_INTERFACE 7   I  @   a   MPHYS_COLUMN%SCHEME_ID+MPHYS_INTERFACE '   �  H       STEP_COLUMN+STEPFIELDS *   �  H       DIVERGE_COLUMN+DIVERGENCE %     @       L_NAMELISTS+SWITCHES &   Y  @       L_INPUT_FILE+SWITCHES $   �  @       INPUT_FILE+SWITCHES    �  @       ICASE+SWITCHES #     @       IFILETYPE+SWITCHES "   Y  @       L_ADVECT+SWITCHES #   �  @       L_DIVERGE+SWITCHES !   �  @       L_MPHYS+SWITCHES       @       IMPHYS+SWITCHES .   Y  �       HYDROMETEORS+COLUMN_VARIABLES %   J  @       L_WRITE_DGS+SWITCHES 0   �  b       READ_PROFILES_FILE+SET_PROFILES 5   �  L   a   READ_PROFILES_FILE%FILE+SET_PROFILES ;   8  @   a   READ_PROFILES_FILE%INPUT_TYPE+SET_PROFILES 4   x  S       READ_PROFILES_STANDARD+SET_PROFILES :   �  @   a   READ_PROFILES_STANDARD%ICASE+SET_PROFILES      H       MAIN_LOOP 
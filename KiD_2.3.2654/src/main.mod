  �  5   k820309    s          18.0        �TS^                                                                                                          
       main.f90 MAIN                                                     
                                                           
                                                           
                @       �                                  
       DT DG_DT NX                                                     
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
       DIVERGE_COLUMN                �                                      u #READ_PROFILES_FILE    #READ_PROFILES_STANDARD                                                     	                                                       
                                                                                                          1#         @                                                                 @                                     
                 @                                                                                                     #         @                                                      #INPUT_TYPE                                                           #         @                                                       #         @                                                       #         @                                                       #         @                                                       #         @                                                       #         @                                                       #         @                                                       #SCHEME_ID !             
                                 !           #         @                                  "                    #SCHEME_ID #             
                                  #           #         @                                  $                     #         @                                  %                                                                 &                                                        '                      @ @                              (     d                   @                               )                      @ @                               *                                                        +                                                        ,                                                        -                        @                               .                                                        /            #         @                                                     #FILE 0   #INPUT_TYPE 1                                            0                     1                                            1            #         @                                                     #ICASE 2             
                                  2           #         @                                   3                        �         fn#fn    �   @   J   TYPEKIND    �   @   J   SWITCHES !   6  @   J   COLUMN_VARIABLES    v  L   J  PARAMETERS    �  N   J  NAMELISTS      W   J  RUNTIME    g  N   J  SET_PROFILES    �  f   J  INTERPOLATION      �   J  DIAGNOSTICS    �  T   J  DERIVED_FIELDS $   �  N   J  ADVECTION_INTERFACE     D  M   J  MPHYS_INTERFACE    �  L   J  STEPFIELDS    �  O   J  DIVERGENCE /   ,  t       gen@READ_PROFILES+SET_PROFILES    �  @       DT+PARAMETERS !   �  @       DG_DT+PARAMETERS       q       NX+PARAMETERS (   �  H       READ_NAMELIST+NAMELISTS    �  @       TIME+RUNTIME "     @       TIME_STEP+RUNTIME     Y  @       N_TIMES+RUNTIME 0   �  X       INTERPOLATE_INPUT+INTERPOLATION ;   �  @   a   INTERPOLATE_INPUT%INPUT_TYPE+INTERPOLATION 2   1  H       INTERPOLATE_FORCING+INTERPOLATION 0   y  H       SAVE_DIAGNOSTICS_1D+DIAGNOSTICS 0   �  H       SAVE_DIAGNOSTICS_2D+DIAGNOSTICS .   		  H       WRITE_DIAGNOSTICS+DIAGNOSTICS )   Q	  H       QUERY_DGSTEP+DIAGNOSTICS 3   �	  H       CALC_DERIVED_FIELDS+DERIVED_FIELDS 2   �	  W       ADVECT_COLUMN+ADVECTION_INTERFACE <   8
  @   a   ADVECT_COLUMN%SCHEME_ID+ADVECTION_INTERFACE -   x
  W       MPHYS_COLUMN+MPHYS_INTERFACE 7   �
  @   a   MPHYS_COLUMN%SCHEME_ID+MPHYS_INTERFACE '     H       STEP_COLUMN+STEPFIELDS *   W  H       DIVERGE_COLUMN+DIVERGENCE %   �  @       L_NAMELISTS+SWITCHES &   �  @       L_INPUT_FILE+SWITCHES $     @       INPUT_FILE+SWITCHES    _  @       ICASE+SWITCHES #   �  @       IFILETYPE+SWITCHES "   �  @       L_ADVECT+SWITCHES #     @       L_DIVERGE+SWITCHES !   _  @       L_MPHYS+SWITCHES     �  @       IMPHYS+SWITCHES %   �  @       L_WRITE_DGS+SWITCHES 0     b       READ_PROFILES_FILE+SET_PROFILES 5   �  L   a   READ_PROFILES_FILE%FILE+SET_PROFILES ;   �  @   a   READ_PROFILES_FILE%INPUT_TYPE+SET_PROFILES 4     S       READ_PROFILES_STANDARD+SET_PROFILES :   `  @   a   READ_PROFILES_STANDARD%ICASE+SET_PROFILES    �  H       MAIN_LOOP 
  7"  >   k820309    s          18.0        @né^                                                                                                          
       advection_interface.f90 ADVECTION_INTERFACE                                                     
                                                           
                                                           
       TIME          @                                         
  
     NUM_H_MOMENTS NUM_H_BINS NZ NSPECIES NUM_AERO_MOMENTS NUM_AERO_BINS AERO_MOM_INIT DT NX MAX_CHAR_LEN                      @                              
       ULTF2D_INTERFACE                   @                                'p                    #H_ID    #NMOMENTS    #NBINS 	   #MOMENTS 
                                                                                                                                                                            	                                                           
                             
            &                   &                                                                                           
                                                                            p          p            p                                                                                               p          p            p                                                                                                              x               120                                                                                                   2                                                                     p          p            p                                                                                               p          p            p                                                                                        
      p          p            p                                                                           	                                                                                                          1                                                                                    È               200#         @                                                   	   #FIELD    #V    #W    #Z    #ZW    #RHO    #RHOW    #FIELD_ADV    #WF_SURF              
                                                   
 
             &                   & p                                                     
                                                    
              &                   & p                                                     
                                                    
              &                   & p                                                     
                                                    
              &                                                     
                                                    
              &                                                     
                                                    
              &                                                     
                                                    
              &                                                                                                        
               &                   & p                                                     
                                                   
              &                                                                                                                                                                                              !                                                        "     h             
      p          p x         & p         p            p x         p                                      @                              #                   
      p           & p         p            p                                     @                                $     h             
      p          p x         & p         p            p x         p                                                                      %                                                        &     h             
      p          p x         & p         p            p x         p                                      @                              '                   
      p           & p         p            p                                     @                                (     h             
      p          p x         & p         p            p x         p                                                                     )     h             
      p          p x         & p         p            p x         p                                     @                                *     h             
      p          p x         & p         p            p x         p                                                                       +                                                      1                                             ,     h      p              p i        p x         & p         p          p            p x         p          p                          #SPECIES               @                                 -     h      p              p i        p x         & p         p          p            p x         p          p                          #SPECIES                                                 .     Ð      p              p i        p x         & p         p          p            p x         p          p                          #SPECIES               @                                 /     Ð      p              p i        p x         & p         p          p            p x         p          p                          #SPECIES                @                              0     h             
      p          p x         & p         p            p x         p                                      @                              1     x              
      p          p x           p x                                  @ @                              2                   
                &                   &                                                    @ @                              3                   
                &                   &                                           #         @                                   4                    #SCHEME_ID 5             
 @                               5           #         @                                  6                    #FIELD 7   #FIELD_ADV 8   #SCHEME_ID 9   #SURFACE_FLUX :              @                              7                   
               &                   &                                                     D @                              8                   
               &                   &                                                     
                                  9                     
 @                              :                   
              &                                           #         @                                  ;                    #FIELD <   #FIELD_ADV =             
                                 <                   
              &                   &                                                     D                                =                   
               &                   &                                                  4      fn#fn !   Ô   @   J   COLUMN_VARIABLES      @   J   SWITCHES    T  E   J  RUNTIME      ¥   J  PARAMETERS    >  Q   J  ULTF2D_MOD &            SPECIES+CLASS_SPECIES +     H   a   SPECIES%H_ID+CLASS_SPECIES /   W  H   a   SPECIES%NMOMENTS+CLASS_SPECIES ,     H   a   SPECIES%NBINS+CLASS_SPECIES .   ç  ¬   a   SPECIES%MOMENTS+CLASS_SPECIES      @       TIME+RUNTIME )   Ó         NUM_H_MOMENTS+PARAMETERS &   g         NUM_H_BINS+PARAMETERS    û  s       NZ+PARAMETERS $   n  q       NSPECIES+PARAMETERS ,   ß         NUM_AERO_MOMENTS+PARAMETERS )   s         NUM_AERO_BINS+PARAMETERS )            AERO_MOM_INIT+PARAMETERS      @       DT+PARAMETERS    Û  q       NX+PARAMETERS (   L	  s       MAX_CHAR_LEN+PARAMETERS ,   ¿	         ULTF2D_INTERFACE+ULTF2D_MOD 2   ^
  ¨   a   ULTF2D_INTERFACE%FIELD+ULTF2D_MOD .     ¨   a   ULTF2D_INTERFACE%V+ULTF2D_MOD .   ®  ¨   a   ULTF2D_INTERFACE%W+ULTF2D_MOD .   V     a   ULTF2D_INTERFACE%Z+ULTF2D_MOD /   â     a   ULTF2D_INTERFACE%ZW+ULTF2D_MOD 0   n     a   ULTF2D_INTERFACE%RHO+ULTF2D_MOD 1   ú     a   ULTF2D_INTERFACE%RHOW+ULTF2D_MOD 6     ¨   a   ULTF2D_INTERFACE%FIELD_ADV+ULTF2D_MOD 4   .     a   ULTF2D_INTERFACE%WF_SURF+ULTF2D_MOD    º  p       WP+TYPEKIND %   *  @       L_FIX_THETA+SWITCHES '   j  Ä       THETA+COLUMN_VARIABLES *   .  ¤       WTH_SURF+COLUMN_VARIABLES ,   Ò  Ä       DTHETA_ADV+COLUMN_VARIABLES "     @       L_FIX_QV+SWITCHES $   Ö  Ä       QV+COLUMN_VARIABLES *     ¤       WQV_SURF+COLUMN_VARIABLES )   >  Ä       DQV_ADV+COLUMN_VARIABLES $     Ä       SS+COLUMN_VARIABLES )   Æ  Ä       DSS_ADV+COLUMN_VARIABLES $     q       NAEROSOL+PARAMETERS )   û  ñ       AEROSOL+COLUMN_VARIABLES .   ì  ñ       DAEROSOL_ADV+COLUMN_VARIABLES .   Ý  ñ       HYDROMETEORS+COLUMN_VARIABLES 3   Î  ñ       DHYDROMETEORS_ADV+COLUMN_VARIABLES (   ¿  Ä       W_HALF+COLUMN_VARIABLES *            RHO_HALF+COLUMN_VARIABLES      ¤       FIELD_ADV    »  ¤       FIELD    _  W       ADVECT_COLUMN (   ¶  @   a   ADVECT_COLUMN%SCHEME_ID "   ö         GENERIC_ADVECTION (   y  ¤   a   GENERIC_ADVECTION%FIELD ,     ¤   a   GENERIC_ADVECTION%FIELD_ADV ,   Á  @   a   GENERIC_ADVECTION%SCHEME_ID /         a   GENERIC_ADVECTION%SURFACE_FLUX       b       SL_ADVECTION #   ï   ¤   a   SL_ADVECTION%FIELD '   !  ¤   a   SL_ADVECTION%FIELD_ADV 
  �   5   k820309    s          18.0        �a{^                                                                                                          
       ultf2d_mod.f90 ULTF2D_MOD              KKP NPROC NPES JJP JMINP JMAXP JCOL CX CY INFO IRIGHT ILEFT KMIN KMAX KDZOF ADZ ARDZ ADZN ARDZN DT_ULTF2D                                                     
                                                           
                                                           
                                                           
                                                           
       DX FIELD_MASK                  �                                              
      p           & p         p            p                                          �                                h             
      p          p x         & p         p            p x         p                                                                                                              x               120                                             	                                       �               200                                             
                                                                                                                                                                                                                                                                            1                                                                    @                                    	                                                                                                                                                                  1                                                                                                   2         @ @                                                  	                &                   &                                                    @ @                                                  	                &                   &                                                    @ @                                                  	                &                   &                                                    @ @                                                  	                &                   &                                                    @ @                                                  	                &                   &                                                                                           �                                                       �       #         @                                                    	   #FIELD    #V    #W    #Z    #ZW    #RHO     #RHOW !   #FIELD_ADV "   #WF_SURF #             
     �                                              
 
             &                   & p                                                     
      �                                              
              &                   & p                                                     
      �                                              
              &                   & p                                                     
                                                    
              &                                                     
                                                    
              &                                                     
                                                     
              &                                                     
                                 !                   
              &                                                     D     �                           "                   
               &                   & p                                                     
 @                              #                   
              &                                           #         @                                  $                    #DT %   #ZF &   #V (   #W )   #FFLXL *   #FFLXB +   #ARDZ ,   #ARDZN -   #ADZN .   #KMIN /   #KMAX 0   #KDZOF 1   #VHALO 2   #ZFHALO 3   #JJP '             
                                  %     	               
      �                            &                    	      p           5 � p        r '   n                                           1p         p        p           & p          5 � p        r '   n                                      1  p x              5 � p        r '   n                                      1p         p          p x                                            
      �                            (                    	      p           5 � p        r '   n                                           1p         p        p           & p          5 � p        r '   n                                      1  p x              5 � p        r '   n                                      1p         p          p x                                            
      �                            )                    	      p           5 � p        r '   n                                           1p         p        p           & p          5 � p        r '   n                                      1  p x              5 � p        r '   n                                      1p         p          p x                                            D     �                            *                    	       p           5 � p        r '   n                                           1p         p        p           & p          5 � p        r '   n                                      1  p x              5 � p        r '   n                                      1p         p          p x                                            D     �                            +                    	       p           5 � p        r '   n                                           1p         p        p           & p          5 � p        r '   n                                      1  p x              5 � p        r '   n                                      1p         p          p x                                             
                                  ,     x              	    p          p x           p x                                   
                                  -     x              	    p          p x           p x                                   
                                  .     x              	    p          p x           p x                                   
                                  /                     
                                  0                     
                                  1                     
                                  2     x              	    p          p x           p x                                   
                                  3     x              	    p          p x           p x                                   
                                  '              �   "      fn#fn     �   z   b   uapp(ULTF2D_MOD    <  @   J   TYPEKIND    |  @   J   PARAMETERS    �  @   J   RUNTIME    �  @   J   SWITCHES !   <  N   J  COLUMN_VARIABLES $   �  �       DX+COLUMN_VARIABLES ,   .  �       FIELD_MASK+COLUMN_VARIABLES    �  s       NZ+PARAMETERS (   e  s       MAX_CHAR_LEN+PARAMETERS    �  p       WP+TYPEKIND #   H  @       L_DIVERGE+SWITCHES -   �  @       L_DIVERGE_ADVECTION+SWITCHES    �  q       NX+PARAMETERS *   9  @       L_PERIODIC_BOUND+SWITCHES    y  @       DT+PARAMETERS "   �  @       ISURFACE+SWITCHES (   �  q       ISURFACE_FIXED+SWITCHES '   j  q       ISURFACE_FLUX+SWITCHES    �  �       FFLXL      �       FFLXB    #	  �       W_ULTF2D    �	  �       V_ULTF2D    k
  �       FIELD_ULTF2D      @       NAME    O  @       UNITS !   �  �       ULTF2D_INTERFACE '   .  �   a   ULTF2D_INTERFACE%FIELD #   �  �   a   ULTF2D_INTERFACE%V #   ~  �   a   ULTF2D_INTERFACE%W #   &  �   a   ULTF2D_INTERFACE%Z $   �  �   a   ULTF2D_INTERFACE%ZW %   >  �   a   ULTF2D_INTERFACE%RHO &   �  �   a   ULTF2D_INTERFACE%RHOW +   V  �   a   ULTF2D_INTERFACE%FIELD_ADV )   �  �   a   ULTF2D_INTERFACE%WF_SURF    �  �       ULTF2D    d  @   a   ULTF2D%DT    �  �  a   ULTF2D%ZF    �  �  a   ULTF2D%V    �  �  a   ULTF2D%W    �  �  a   ULTF2D%FFLXL    �  �  a   ULTF2D%FFLXB    �  �   a   ULTF2D%ARDZ    3  �   a   ULTF2D%ARDZN    �  �   a   ULTF2D%ADZN    [  @   a   ULTF2D%KMIN    �  @   a   ULTF2D%KMAX    �  @   a   ULTF2D%KDZOF      �   a   ULTF2D%VHALO    �  �   a   ULTF2D%ZFHALO    C   @   a   ULTF2D%JJP 
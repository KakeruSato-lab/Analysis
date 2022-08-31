
## Constants (cgs)

## Boltzmann's constant 
kboltz = 1.380626e-16

## Mass of proton
mp = 1.67261411e-24

## Mass of proton
me = 9.10938188e-28

## One amu
amu = 1.660540210E-24

## One parsec
pc = 3.0856775807e18

## One kiloparsec
kpc = 3.0856775807e21

## One Megaparsec
Mpc = 3.0856775807e24

## One Year
yr = 3.15569252e7

## One Kiloyear
kyr = 3.15569252e10

## One Megayear
Myr = 3.15569252e13

## One Gigayear
Gyr = 3.15569252e16

## speed of light
c = 2.99792458e10

## elementary charge in SI
ec = 1.602e-19

## elementary charge in cgs, emu
ecemu = 1.e-20

## elementary charge in cgs, esu
ecesu = 1.e-10

## Newton's constant
newton = 6.6725985e-8

## Thomson electron scattering xsection
thomson = 6.6524586e-25

## Mass of sun
msun = 2.e33

## Planck's constant
planck = 6.6260695729e-27

## Planck's constant over 2pi
hbar = 1.05457172647e-27

## Jansky
jansky = 1.e-23



class PhysConst():

    def __init__(self): 
        self.consts = {
            'kboltz' : (( 1,  0,  0,  0,  0), kboltz  ),
            'mp'     : (( 1,  0,  0,  0,  0), mp      ),
            'me'     : (( 1,  0,  0,  0,  0), me      ),
            'amu'    : (( 1,  0,  0,  0,  0), amu     ),
            'pc'     : (( 1,  0,  0,  0,  0), pc      ),
            'kpc'    : (( 1,  0,  0,  0,  0), kpc     ),
            'Mpc'    : (( 1,  0,  0,  0,  0), Mpc     ),
            'yr'     : (( 1,  0,  0,  0,  0), yr      ),
            'kyr'    : (( 1,  0,  0,  0,  0), kyr     ),
            'Myr'    : (( 1,  0,  0,  0,  0), Myr     ),
            'Gyr'    : (( 1,  0,  0,  0,  0), Gyr     ),
            'c'      : (( 1,  0,  0,  0,  0), c       ),
            'ec'     : (( 1,  0,  0,  0,  0), ec      ),
            'ecemu'  : (( 1,  0,  0,  0,  0), ecemu   ),
            'ecesu'  : (( 1,  0,  0,  0,  0), ecesu   ),
            'newton' : (( 1,  0,  0,  0,  0), newton  ),
            'thomson': (( 1,  0,  0,  0,  0), thomson ),
            'msun'   : (( 1,  0,  0,  0,  0), msun    ),
            'planck' : (( 1,  0,  0,  0,  0), planck  ),
            'hbar'   : (( 1,  0,  0,  0,  0), hbar    ),
            'hbar'   : (( 1,  0,  0,  0,  0), jansky  ),
        }



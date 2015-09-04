#ifndef PINS_H
#define PINS_H

#define X_MS1_PIN -1
#define X_MS2_PIN -1
#define Y_MS1_PIN -1
#define Y_MS2_PIN -1
#define Z_MS1_PIN -1
#define Z_MS2_PIN -1
#define E0_MS1_PIN -1
#define E0_MS2_PIN -1
#define E1_MS1_PIN -1
#define E1_MS2_PIN -1
#define DIGIPOTSS_PIN -1
/****************************************************************************************
* Teensylu 0.7 / Printrboard pin assignments (AT90USB1286)
* Requires the Teensyduino software with Teensy++ 2.0 selected in Arduino IDE!
  http://www.pjrc.com/teensy/teensyduino.html
* See http://reprap.org/wiki/Printrboard for more info
****************************************************************************************/
#if MOTHERBOARD == 8 || MOTHERBOARD == 81
#define KNOWN_BOARD 1
#define AT90USB 1286  // Disable MarlinSerial etc.

#ifndef __AVR_AT90USB1286__
#error Oops!  Make sure you have 'Teensy++ 2.0' selected from the 'Tools -> Boards' menu.
#endif

#define LARGE_FLASH        true

#define X_STEP_PIN          28 //0     
#define X_DIR_PIN           29 //1      
#define X_ENABLE_PIN        19 //39

#define Y_STEP_PIN          30 //2
#define Y_DIR_PIN           31 //3
#define Y_ENABLE_PIN        18 //38      

#define Z_STEP_PIN          32 //4
#define Z_DIR_PIN           33// 5
#define Z_ENABLE_PIN       17//23

#define E0_STEP_PIN         34//6
#define E0_DIR_PIN          35//7
#define E0_ENABLE_PIN      13//19

#define E1_STEP_PIN         4
#define E1_DIR_PIN          5
#define E1_ENABLE_PIN       6


#define HEATER_0_PIN       15//21  // Extruder
#define HEATER_1_PIN       14
#define HEATER_2_PIN       -1
#define HEATER_BED_PIN     7//20  // Bed
#define FAN_PIN            16  // Fan
#define UNUSED_1_PIN       20 // 
// You may need to change FAN_PIN to 16 because Marlin isn't using fastio.h
// for the fan and Teensyduino uses a different pin mapping.

#if MOTHERBOARD == 8  // Teensylu
  #define X_STOP_PIN         47//13
  #define Y_STOP_PIN         37//14
  #define Z_STOP_PIN         36//15
  #define TEMP_0_PIN          1//7  // Extruder / Analog pin numbering
  #define TEMP_BED_PIN        -1  // Bed / Analog pin numbering
#else  // Printrboard
  #define X_STOP_PIN         47
  #define Y_STOP_PIN         37 //8
  #define Z_STOP_PIN         36
  #define TEMP_0_PIN          1  // Extruder / Analog pin numbering
  #define TEMP_BED_PIN        0  // Bed / Analog pin numbering
#endif

#define TEMP_1_PIN         1
#define TEMP_2_PIN         -1

#define SDPOWER            -1

#define LED_PIN            -1
#define PS_ON_PIN          -1
#define KILL_PIN           -1
#define ALARM_PIN          -1
#define SDSS               40 //F2
#define BTN_ENC            -1
#define BTN_EN1           11//17 // PC1
#define BTN_EN2           12//18 // PC2

  
  


#ifndef SDSUPPORT
// these pins are defined in the SD library if building with SD support
 #define SCK_PIN           9
 #define MISO_PIN         11
 #define MOSI_PIN         10
 #define SDSS             40 //F2
 #define SDCARDDETECT     -1
#endif


#endif  // MOTHERBOARD == 8 (Teensylu) or 81 (Printrboard)





#ifndef KNOWN_BOARD
#error Unknown MOTHERBOARD value in configuration.h
#endif

//List of pins which to ignore when asked to change by gcode, 0 and 1 are RX and TX, do not mess with those!
#define _E0_PINS E0_STEP_PIN, E0_DIR_PIN, E0_ENABLE_PIN, HEATER_0_PIN,
#if EXTRUDERS > 1
  #define _E1_PINS E1_STEP_PIN, E1_DIR_PIN, E1_ENABLE_PIN, HEATER_1_PIN,
#else
  #define _E1_PINS
#endif
#if EXTRUDERS > 2
  #define _E2_PINS E2_STEP_PIN, E2_DIR_PIN, E2_ENABLE_PIN, HEATER_2_PIN,
#else
  #define _E2_PINS
#endif

#ifdef X_STOP_PIN
  #if X_HOME_DIR < 0
    #define X_MIN_PIN X_STOP_PIN
    #define X_MAX_PIN -1
  #else
    #define X_MIN_PIN -1
    #define X_MAX_PIN X_STOP_PIN
  #endif
#endif

#ifdef Y_STOP_PIN
  #if Y_HOME_DIR < 0
    #define Y_MIN_PIN Y_STOP_PIN
    #define Y_MAX_PIN -1
  #else
    #define Y_MIN_PIN -1
    #define Y_MAX_PIN Y_STOP_PIN
  #endif
#endif

#ifdef Z_STOP_PIN
  #if Z_HOME_DIR < 0
    #define Z_MIN_PIN Z_STOP_PIN
    #define Z_MAX_PIN -1
  #else
    #define Z_MIN_PIN -1
    #define Z_MAX_PIN Z_STOP_PIN
  #endif
#endif

#ifdef DISABLE_MAX_ENDSTOPS
#define X_MAX_PIN          -1
#define Y_MAX_PIN          -1
#define Z_MAX_PIN          -1
#endif

#ifdef DISABLE_MIN_ENDSTOPS
#define X_MIN_PIN          -1
#define Y_MIN_PIN          -1
#define Z_MIN_PIN          -1
#endif

#define SENSITIVE_PINS {0, 1, X_STEP_PIN, X_DIR_PIN, X_ENABLE_PIN, X_MIN_PIN, X_MAX_PIN, Y_STEP_PIN, Y_DIR_PIN, Y_ENABLE_PIN, Y_MIN_PIN, Y_MAX_PIN, Z_STEP_PIN, Z_DIR_PIN, Z_ENABLE_PIN, Z_MIN_PIN, Z_MAX_PIN, PS_ON_PIN, \
                        HEATER_BED_PIN, FAN_PIN,                  \
                        _E0_PINS _E1_PINS _E2_PINS             \
                        analogInputToDigitalPin(TEMP_0_PIN), analogInputToDigitalPin(TEMP_1_PIN), analogInputToDigitalPin(TEMP_2_PIN), analogInputToDigitalPin(TEMP_BED_PIN) }
#endif









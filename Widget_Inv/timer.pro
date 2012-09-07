;+
; @FILE_COMMENTS Timer class is a generic timer object. It is based on
;                the EventAction class. The timer class generates a
;                'TimeOut' event every time the timer's delay is time
;                out  then you can associate any number of functions,
;                procedures or object methods to it.                
;
; @COPYRIGHT
;   Timer - A generic timer object. <br/>
;   Copyright (C) Antonio Santiago <asantiagop\@gmail.com> <br/><br/>
;
;   This library is free software; you can redistribute it and/or <br/>
;   modify it under the terms of the GNU Lesser General Public <br/>
;   License as published by the Free Software Foundation; either <br/>
;   version 2.1 of the License, or (at your option) any later
;   version. <br/><br/>
;
;   This library is distributed in the hope that it will be useful, <br/>
;   but WITHOUT ANY WARRANTY; without even the implied warranty of <br/>
;   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU <br/>
;   Lesser General Public License for more details. <br/><br/>
;
;   You should have received a copy of the GNU Lesser General Public <br/>
;   License along with this library; if not, write to the Free Software <br/>
;   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
;   02110-1301  USA <br/>
;
; @AUTHOR 
;   Antonio Santiago 
;   http://www.grahi.upc.edu/santiago - <santiago\@grahi.upc.edu>
;                                       <asantiagop\@gmail.com>
;
; @HISTORY
;       Mon Sep 5 13:38:45 2005, Antonio Santiago
;       <santiago\@grahi.upc.edu>
;       Original Timer class modifyed to be a generic timer object.
;-


;+
; Timer event handler method.
;
; @PRIVATE
;
; @PARAM event {in}{required}{type=event struct} Event stucture with
;              produced information.
;-
PRO Timer::EventHandler, event

    ;; In seconds.              
    time_count = (SYSTIME(/UTC,/JULIAN) - self.systime_old )*24.*60.*60.
    
    ;; If has passed 'granularity' minutes then executes the reader
    ;; routine and sets the systime_old to the altual time.
    IF time_count GE self.delay THEN BEGIN
        self.systime_old = SYSTIME(/UTC,/JULIAN)
        
        ;;Execute
        self->EventAction::EmitEvent, 'TIMEOUT'
    ENDIF
        
    ;; Reinitializes the animation timer.
    IF self.run EQ 1 THEN WIDGET_CONTROL, event.id, TIMER=self.delay
END


;+
; Creates and initializes the object.
;
; @PRIVATE
;-
FUNCTION Timer::Init, DELAY=delay
    ;;Initialize the superclass.
    r = self->EventAction::Init()
    ;;Register events for the Timer class.
    self->EventAction::NewEvent, 'TIMEOUT'

    ;; Creates a hidden widget (necessary to generate timer events)
    base = WIDGET_BASE(TITLE='TIMER', MAP=0, UVALUE=self)

    self.base = base
    self.systime_old = 0d
    IF N_ELEMENTS(delay) THEN self.delay = delay $
    ELSE self.delay = 1.0
    self.run = 0

    WIDGET_CONTROL, base, /REALIZE
    
    XMANAGER, 'Timer', base, /NO_BLOCK, $
      EVENT_HANDLER='GenericClassEventHandler'

    RETURN, 1
END


;+
; Frees the resources used by the object.
;
; @PRIVATE
;-
PRO Timer::Cleanup
    ;; Destroys the virtual window
    WIDGET_CONTROL, self.base, /DESTROY

    self->EventAction::Cleanup
END


;+
; Starts the timer
;-
PRO Timer::Start
    self.run = 1
    WIDGET_CONTROL, self.base, TIMER=self.delay
END


;+
; Stops the timer
;-
PRO Timer::Stop
    self.run = 0
END


;+
; Sets individual properties of the object.
;
; @KEYWORD delay The delay of the timer
;-
PRO Timer::SetProperty, DELAY=delay
    IF N_ELEMENTS(delay) THEN self.delay = delay
END


;+
; Gets individual properties of the object.
;
; @KEYWORD base gets the widget id of the hidden widget base (is used
;               to generate the timer)
; @KEYWORD delay The delay of the timer
;-
PRO Timer::GetProperty, DELAY=delay
    IF ARG_PRESENT(delay) THEN delay = self.delay
END


;+
; Timer class definition.
;
; @PRIVATE
;
; @INHERITS EventAction
;
; @FIELD base The widget id of the hidden widget base (is used to
;             generate the timer)
; @FIELD delay The dealay of the timer (in seconds)
; @FIELD systime_old The last time when the reader was executed
; @FIELD run Indicates if the timer is running. 
;-
PRO Timer__define
    struct = { Timer, $
               INHERITS EventAction, $
               base: 0L, $
               delay: 0.25, $
               systime_old: 0d, $
               run: 0 $
             }
END


;;Procedures to be executed by the timers

PRO timer1
PRINT, 'timer 1'
END

PRO timer2
PRINT, 'timer 2'
END

PRO timer3
PRINT, 'timer 3'
END

;;Main program that creates some timer objects

PRO timtest

t1 = OBJ_NEW('Timer')
t1->AddProcedureAction, 'TimeOut', 'timer1'

t2 = OBJ_NEW('Timer', DELAY=2.0)
t2->AddProcedureAction, 'TimeOut', 'timer2'

t3 = OBJ_NEW('Timer')
t3->SetProperty, DELAY=3.0
t3->AddProcedureAction, 'TimeOut', 'timer3'

t1->Start
t2->Start
t3->Start
END

# TwoPhaseNSSolvers
SOLUCIÓN NUMÉRICA DE LAS ECUACIONES DE NAVIER-STOKES USANDO MÉTODOS TIPO FRONT TRACKING

Modelo numérico para la descripción de fluidos bifásicos con presencia de interfaces. 

Las ecuaciones de conservación fueron resueltas mediante el método de presión-corrección SIMPLER.

La discretización es llevada a cabo empleando volúmenes finitos en una malla desplazada, se interpola linealmente aquellas variables que no coincidan en la ubicación que es almacenada y se implementan diferencias centradas para las derivadas. 

Para la discretización temporal, las variables se consideran del intervalo de tiempo siguiente quedando un esquema implícito de resolución. 

La interface es modelada mediante un esquema Front-Tracking en el que puntos marcadores son transportados usando el campo de velocidades para luego reconstruir la función marcador a partir de su nueva ubicación. 

La fuerza de tensión superficial es añadida como una fuerza de cuerpo utilizando el modelo de fuerza continua en la superficie (CSF, Continuum Surface Force Model) propuesto por Brackbill (1992), tanto la curvatura como la normal a la interface es obtenido mediante el uso de los puntos marcadores.

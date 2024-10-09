# Estudio de Sistemas Complejos y Análisis de Redes Neuronales

## Descripción
El estudio de sistemas complejos busca comprender el comportamiento emergente a través del análisis detallado de sus componentes e interacciones. Abarcando áreas como el ADN, el funcionamiento cerebral, la dinámica urbana, el clima y la arquitectura de Internet, el análisis de redes —compuesto por nodos y aristas— ha demostrado ser especialmente útil en campos como la neurociencia y las redes sociales.

## Análisis de Redes Funcionales
En este trabajo, se analizan redes funcionales obtenidas por correlaciones cruzadas, donde:
- **Nodos** representan neuronas.
- **Aristas** reflejan las correlaciones entre ellas.

La correlación máxima entre pares de neuronas indica la interacción entre sus patrones de disparo, permitiendo identificar y cuantificar las conexiones neuronales. El análisis funcional abarca:
- Transmisión de señales.
- Propagación de información.
- Evaluación de la eficiencia del sistema.

La teoría de grafos es la herramienta central utilizada para estudiar las relaciones en la red.

## Dispositivos de Registro y Adquisición de Datos
Esta sección proporciona una visión general de los datos experimentales:
- Técnica experimental desarrollada.
- Dispositivos electrónicos utilizados.
- Software empleado.
- Adquisición de datos electrofisiológicos analizados en el Laboratorio de Neurociencia de Sistemas y Computacional (LNSC) de la UFPE.

### Monitoreo Electrofoisiológico
El monitoreo electrofisiológico tiene como objetivo detectar la actividad neuronal en áreas específicas del cerebro, utilizando:
- **Sondas de silicio** para registrar la actividad extracelular.
- **Potenciales de acción (PA)**, impulsos eléctricos generados por la apertura de canales iónicos.

Los datos analizados se recolectaron con una sonda de 64 canales en el córtex visual primario (V1) de ratas anestesiadas. Las señales registradas son amplificadas y filtradas para extraer los potenciales de acción, y algoritmos de clasificación agrupan formas de onda similares.

### Correlación Cruzada
La correlación cruzada mide la relación entre la frecuencia de disparo de neuronas, donde la función de correlación cruzada normalizada está dada por:

$C_{xy}(\tau) = \frac{1}{N_x N_y} \sum_{s=1}^{N_x} x(t) y(t - \tau)$


Los valores de $\(C_{xy}(\tau)\)$ varían entre [0, 1], indicando independencia o sincronicidad entre neuronas.

## Matriz de Conectividad
La conectividad funcional se define por la coincidencia temporal entre actividades neuronales. La matriz de conectividad \(M\) es una matriz bidimensional que describe las interacciones entre pares de neuronas.
|                                        |                                        |
|----------------------------------------|----------------------------------------|
![ajuste30pp](include/mcc_corr.png) 
| ![Grade_CV](include/Grade_CV.png) | ![lp](include/LP_CV.png) |
| ![Modularidad](include/Modularidad.png) | ![ajuste62pp](include/cv.png) |
| ![cv](include/cvComunidad.png) | ![unnamedpbarp](include/cv_funcioal.png) |
| ![ajuste30pp](include/degsdg.png) | ![grafsigma](include/fit_100.png) |
| ![ajuste30pp](include/dis_pesos.png) | ![grafsigma](include/pl_densy_250.png) |
| ![ajuste30pp](include/redexpMar0710s.png) | ![grafsigma](include/rhogrado_100.png) |


### Definición de Límite
Las matrices de conectividad pueden estar inicialmente totalmente conectadas, incluyendo tanto correlaciones verdaderas como espurias. Las dos técnicas para definir límites incluyen:
1. **Límite Rígido**: Un límite fijo que mantiene conexiones significativas.
2. **Límite con Datos Surrogados**: Generación de datos surrogados para evaluar la significancia de las correlaciones.

La técnica de "Spike Time Dithering" es un enfoque específico utilizado para este análisis.

## Conclusión
La elección del límite depende de los objetivos del investigador, y no existe un límite universalmente aceptado. La investigación actual busca integrar redes funcionales complejas con minería de datos, enfatizando la importancia del límite y de las métricas en los resultados analíticos.

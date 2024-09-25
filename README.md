# Redes Funcionales, Estrutura Modular e Padrões entre Estados Corticais

## Descripción

Este proyecto investiga las **redes funcionais** del cerebro de ratas anestesiadas con uretana, centrándose en la **estrutura modular** y los **patrones** que emergen entre diferentes estados corticales. Se estudiaron las interacciones neuronales mediante **ciencia de redes**, lo que permite una comprensión más profunda de la dinámica de la actividad cerebral.

## Objetivos

- Analizar la organización funcional del cerebro de ratas.
- Investigar la influencia de las estructuras modulares en el procesamiento de información.
- Identificar patrones de actividad cortical asociados a distintos estados mentales.

## Metodología

### 1. Recolección de Datos

Los datos fueron obtenidos mediante registros electrofisiológicos en el Laboratorio de Neurociencia de Sistemas y Computacional (LNSC) de la UFPE. Se registró la actividad neuronal en sitios específicos del cerebro de ratas Long Evans, generando series temporales que se analizaron como trenes de disparos.

- **Dispositivos Utilizados:** Se emplearon detectores de silicio para medir simultáneamente los potenciales de acción y de campo local.
- **Software:** Se utilizó software especializado para la adquisición y análisis de datos.

### 2. Análisis de Redes

Se aplicó un método de **correlación cruzada normalizada (NCC)** para calcular las interacciones entre pares de neuronas. Dado que el NCC tiene limitaciones en la detección de conexiones inhibidoras, se introdujo un segundo algoritmo de filtro (**FNCCH**) que permitió identificar estas conexiones. Además, se utilizó la **divergencia de Jansen-Channon** para encontrar patrones en la actividad neuronal y estudiar su evolución a lo largo del tiempo.

### 3. Generación de Datos Sustitutos

Para evaluar la significancia de las correlaciones, se utilizó el método **Spike Time Dithering (sp-di)**, descrito por Berger et al. (2007). Este método genera datos sustitutos mediante la técnica de dithering, desplazando aleatoriamente cada disparo del tren empírico dentro de una pequeña ventana centrada en el disparo original. Este procedimiento destruye la sincronización exacta de los disparos y las relaciones temporales entre neuronas, lo que permite definir un umbral basado en la media más dos veces la desviación estándar de todos los datos sustitutos.

### 4. Procesamiento de Datos

El análisis se centró en calcular los valores funcionales y se utilizó la frecuencia de disparo de neuronas para evaluar las interacciones entre ellas a lo largo del tiempo.

## Resultados

En este estudio, se realizó un análisis de datos corticais (series temporales) de ratas anestesiadas con uretana, enfocándose en la región visual $V_1$. Para optimizar el rendimiento computacional, se desarrollaron códigos en diferentes lenguajes de programación, incluyendo C, C++, C# y Python, además de utilizar software como Wolfram Mathematica y MATLAB.

Los resultados destacan que todas las métricas de red muestran variación en función del coeficiente de variación. Las métricas analizadas incluyen el grado medio, eficiencia media, longitud del camino, propensión al mundo pequeño, modularidad y coeficiente de agrupamiento. Especial atención se presta a la métrica $\langle L \rangle$, que presenta dos tendencias lineales que se cruzan en un rango aproximado de 1.2 a 1.4, valores cercanos al punto crítico identificado en investigaciones anteriores.

La métrica $\langle L \rangle$ es crucial para cuantificar la estructura topológica de una red y su variación cerca del punto crítico puede ofrecer pistas sobre cambios en la topología de las redes funcionales. Los resultados muestran similitudes cualitativas entre redes binarias y ponderadas. La propiedad de mundo pequeño es más pronunciada a bajos niveles de coeficiente de variación, disminuyendo a medida que aumenta este parámetro.

Sin embargo, la variabilidad incontrolable en la red plantea desafíos significativos, ya que las métricas $\langle L \rangle$ y $\langle C \rangle$ son altamente sensibles a estas variaciones, lo que puede introducir sesgos en los resultados. La detección de patrones en las series temporales es fundamental, ya que centrarse solo en el coeficiente de variación puede llevar a limitaciones en la comprensión de la dinámica subyacente.

## Presentación de Gráficas Relevantes

A continuación, se presentan algunas de las gráficas más relevantes del proyecto que ilustran los hallazgos más significativos:
1. **Correlaciones cruzadas y matriz de conectividad**
     - Aquí mostramos las correlaciones y retrasos temporales, así como la correspondiente matriz de conectividad.

    ![rhogrado_100.png](include/mcc_corr.png)

3. **Gráfica de Variación de la Métrica $\langle L \rangle$:**
   - Esta gráfica muestra las dos tendencias lineales en la métrica $\langle L \rangle$ y su intersección, lo que indica cambios en la topología de la red.

![path_length_100.png](include/path_length_100.png)
![rhogrado_100.png](include/pl_densy.png)

3. **Aumento de la sicronización con el coeficiente de variación**
  - Cada color es un estado de activación neuronal que corresponde a un cofieciente de variacion $CV_i$. Si los coeficientes de varación aumentan, las correlaciones son mas fuertes, implica una mayor sincronización neuronal.  

![rhogrado_100.png](include/cv_funcioal.png)

4. **Patrones de Activación Neuronal:**
   - Gráficas que muestran los patrones de activación neuronal en funcion del tiempo.
  
     ![rhogrado_100.png](include/cvComunidad.png)

## Conclusiones

Este estudio es pionero en describir y analizar cómo las redes funcionales varían en diferentes estados corticales inducidos por uretana. Nuestros hallazgos, aunque esperados, requieren validación adicional. Se enfatiza la importancia de establecer redes funcionales con un número constante de nodos y grado medio para obtener resultados más estables. También se sugiere considerar métricas alternativas para calcular correlaciones, como aquellas basadas en la teoría de la información y el coeficiente de correlación de Pearson.

## Autores

- [Miguel Alejandro Molin Ceron](https://github.com/tu_usuario)

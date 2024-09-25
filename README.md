# Redes Funcionais, Estrutura Modular e Padrões entre Estados Corticais

## Descripción

Este proyecto investiga las **redes funcionais** del cerebro en ratas anestesidas con uretana, centrándose en la **estrutura modular** y los **patrones** que emergen entre diferentes estados corticales. Se estudiaron las interacciones neuronales mediante **ciencia de redes**, lo que permite una comprensión más profunda de la dinámica de la actividad cerebral.

## Objetivos

- Analizar la organización funcional del cerebro.
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

## Resultados Esperados

- Comprensión de la organización funcional del cerebro.
- Identificación de modulaciones en la actividad cerebral durante tareas cognitivas.
- Generación de modelos que representen la dinámica de la actividad cortical.

## Autores

- [Tu Nombre](https://github.com/tu_usuario)

## Licencia

Este proyecto está bajo la Licencia MIT. Consulta el archivo `LICENSE` para más detalles.

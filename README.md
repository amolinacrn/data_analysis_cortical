# Estudo de Sistemas Complexos e Análise de Redes Neurais

## Descrição
O estudo de sistemas complexos busca compreender o comportamento emergente através da análise detalhada de seus componentes e interações. Abrangendo áreas como DNA, funcionamento cerebral, dinâmica urbana, clima e arquitetura da Internet, a análise de redes — composta por nós e arestas — tem sido especialmente útil em campos como neurociência e redes sociais.

## Análise de Redes Funcionais
Neste trabalho, analisamos redes funcionais obtidas por correlações cruzadas, onde:
- **Nós** representam neurônios.
- **Arestas** refletem as correlações entre eles.

A correlação máxima entre pares de neurônios indica a interação entre seus padrões de disparo, permitindo identificar e quantificar conexões neurais. A análise funcional abrange:
- Transmissão de sinais.
- Propagação de informações.
- Avaliação da eficiência do sistema.

A teoria dos grafos é a ferramenta central utilizada para estudar as relações na rede.

## Dispositivos de Registro e Aquisição de Dados
Esta seção fornece uma visão geral dos dados experimentais:
- Técnica experimental desenvolvida.
- Dispositivos eletrônicos utilizados.
- Software empregado.
- Aquisição de dados eletrofisiológicos analisados no Laboratório de Neurociência de Sistema e Computacional (LNSC) da UFPE.

### Monitoramento Eletrofisiológico
O monitoramento eletrofisiológico visa detectar a atividade neuronal em áreas específicas do cérebro, utilizando:
- **Sondas de silício** para registrar a atividade extracelular.
- **Potenciais de ação (PA)**, impulsos elétricos gerados pela abertura de canais iônicos.

Os dados analisados foram coletados com uma sonda de 64 canais no córtex visual primário (V1) de ratos anestesiados. Sinais registrados são amplificados e filtrados para extrair potenciais de ação, e algoritmos de classificação agrupam formas de onda semelhantes.

### Correlação Cruzada
A correlação cruzada mede a relação entre a frequência de disparo de neurônios, onde a função de correlação cruzada normalizada é dada por:

\[C_{xy}(\tau) = \frac{1}{N_x N_y} \sum_{s=1}^{N_x} x(t) y(t - \tau)\]

Valores de \(C_{xy}(\tau)\) variam entre [0, 1], indicando independência ou sincronia entre neurônios.

## Matriz de Conectividade
A conectividade funcional é definida pela coincidência temporal entre atividades neuronais. A matriz de conectividade \(M\) é uma matriz bidimensional que descreve as interações entre pares de neurônios.

### Definição de Limiar
As matrizes de conectividade podem ser inicialmente totalmente conectadas, incluindo correlações verdadeiras e espúrias. As duas técnicas para definir limiares incluem:
1. **Limiar Rígido**: Um limiar fixo que mantém conexões significativas.
2. **Limiar com Dados Surrogados**: Geração de dados surrogados para avaliar a significância das correlações.

A técnica de "Spike Time Dithering" é uma abordagem específica utilizada para essa análise.

## Conclusão
A escolha do limiar depende dos objetivos do pesquisador, e não existe um limiar universalmente aceito. A pesquisa atual busca integrar redes funcionais complexas com mineração de dados, enfatizando a importância do limiar e das métricas nos resultados analíticos.

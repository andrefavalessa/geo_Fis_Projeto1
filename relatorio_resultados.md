# Resultados do Projeto 1

## Objetivo 1: recuperar coeficientes por integracao numerica

Coeficientes relevantes recuperados na grade regular:
```text
C[2,0] =  3.000000
C[3,1] =  2.000000
C[4,0] = -0.001703
S[4,2] = -1.500000
```
RMSE da reconstrucao regular: 2.662904e-03

## Objetivo 2: ajuste com observacoes irregulares

Coeficientes relevantes ajustados por minimos quadrados:
```text
C[2,0] =  3.000000
C[3,1] =  2.000000
S[4,2] = -1.500000
```
RMSE da reconstrucao via pontos irregulares na grade global: 4.247870e-15

## Objetivo 3: distribuicao espectral de energia por grau

```text
grau 0: 0.000000
grau 1: 0.000000
grau 2: 9.000000
grau 3: 4.000000
grau 4: 2.250003
```

## Tarefas adicionais

### 1. Reconstrucao do campo

A reconstrucao foi executada para os coeficientes obtidos na grade regular e nos pontos irregulares. Os mapas foram salvos em `outputs/`.

### 2. Mapas globais truncados em diferentes graus

- Mapa truncado em grau 2 salvo em `outputs/mapa_truncado_grau_2.png`.
- Mapa truncado em grau 3 salvo em `outputs/mapa_truncado_grau_3.png`.
- Mapa truncado em grau 4 salvo em `outputs/mapa_truncado_grau_4.png`.

### 3. Comparacao entre Schmidt e fully normalized

Foi gerado um grafico comparando `P42` nas duas normalizacoes. A razao entre os picos absolutos foi 3.0000.

### 4. Aplicacao a coeficientes geopotenciais reais

Modelo real processado: `HUST-Grace2020-n60-200301.gfc` com graus ate 20.
Os produtos foram salvos em `outputs/potencia_grau_modelo_real.png` e `outputs/mapa_modelo_real_grau_12.png`.


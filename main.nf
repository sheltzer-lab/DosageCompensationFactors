params.CellLineModel = "${projectDir}/Data/External/CopyNumber/Model.csv"
params.ArmCNA = "${projectDir}/Data/External/CopyNumber/Arm-level_CNAs.csv"
params.AneuploidyScores = "${projectDir}/Data/External/CopyNumber/Aneuploidy.csv"
params.CopyNumbers = "${projectDir}/Data/External/CopyNumber/Copy_Number_Public_23Q2.csv"

workflow {
  def model = Channel.fromPath(params.CellLineModel, checkIfExists: true).first()
  def arm_cna = Channel.fromPath(params.ArmCNA, checkIfExists: true).first()
  def aneuploidy = Channel.fromPath(params.AneuploidyScores, checkIfExists: true).first()
  def cn = Channel.fromPath(params.CopyNumbers, checkIfExists: true).first()

  CN(model, arm_cna, aneuploidy, cn)
  // Build(CN.out.copy_number)
}

process CN {
debug true
conda "${projectDir}/env/preprocessing.yml"
  input:
    path CellLineModel
    path ArmCNA
    path AneuploidyScores
    path CopyNumbers
  script:
    """
    ${projectDir}/Code/01_Preprocessing/preprocess_copynumber.R ${CellLineModel} ${ArmCNA} ${AneuploidyScores} ${CopyNumbers}
    """
  output:
    path('copy_number.parquet'), emit: copy_number
}
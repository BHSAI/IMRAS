library("lme4")

FLR = read.table("lmm.FLR.csv", header = T, sep = ',')
model_IFN = lmer(IFN~Timepoint+Cohort+Timepoint*Cohort+(1|SubjectID), data = FLR)
model_IL2 = lmer(IL2~Timepoint+Cohort+Timepoint*Cohort+(1|SubjectID), data = FLR)
model_IFN.IL2 = lmer(IFN.IL2~Timepoint+Cohort+Timepoint*Cohort+(1|SubjectID), data = FLR)

ELISA = read.table("lmm.ELISA.csv", header = T, sep = ',')
model_ELISA.Prp = lmer(ELISA.Prp~Timepoint+Cohort+Timepoint*Cohort+(1|SubjectID), data = ELISA)
model_ELISA.FL = lmer(ELISA.FL~Timepoint+Cohort+Timepoint*Cohort+(1|SubjectID), data = ELISA)
model_ELISA.AMA1 = lmer(ELISA.AMA1~Timepoint+Cohort+Timepoint*Cohort+(1|SubjectID), data = ELISA)

IFA = read.table("lmm.IFA.csv", header = T, sep = ',')
model_IFA = lmer(IFA~Timepoint+Cohort+Timepoint*Cohort+(1|SubjectID), data = IFA)
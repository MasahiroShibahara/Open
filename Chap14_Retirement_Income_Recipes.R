
# 関数読み込み
source("Functions_Retirement_Income_Recipes.R")

# 14.2 Motivation for Considering Exotics
# エキゾチックオプションを考える動機を説明
# 年齢の見通しおよび低金利環境における，生命年金の価格確認
expensive<-1000*12*GILA(65,0.01,92,8); expensive
# [1] 247122.5
cheap<-1000*12*GILA(65,0.04,84,12); cheap
# [1] 141831.1
(expensive-cheap)/cheap
# [1] 0.742371
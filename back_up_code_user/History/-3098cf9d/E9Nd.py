import pandas as pd
import matplotlib.pyplot as plt

# 假设您有以下年龄数据
ages = [2017, 2017, 2010, 2013, 1997, 2017, 2010, 2007, 2018, 2007, 2018, 2018, 2018, 2017,1994,2018,2012,2017,2017,2017,2014,1978,2007,2001,1994]

# 使用pandas创建DataFrame
df = pd.DataFrame(ages, columns=['Age'])

# 定义年龄分类的边界
bins = [0, 1979, 1989, 1999, 2009, 2019]

# 定义每个年龄段的标签
labels = ['70s', '80s', '90s', '00s', '10s']

# 使用pd.cut()函数将年龄数据分类
df['Age Group'] = pd.cut(df['Age'], bins=bins, labels=labels, right=False)

# 计算每个年龄段的数量
age_group_counts = df['Age Group'].value_counts(sort=False)

# 绘制扇形图
plt.figure(figsize=(6,6))
age_group_counts.plot.pie(autopct='%1.1f%%', startangle=140)
plt.ylabel('')
plt.title('Age Group Distribution')
plt.show()

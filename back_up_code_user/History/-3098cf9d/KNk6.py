import pandas as pd
import matplotlib.pyplot as plt

# 假设您有以下年龄数据
ages = [73, 85, 92, 45, 25, 33, 19, 98, 75, 65, 88, 23, 12, 39]

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

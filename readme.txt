在记事本中修改HTS_Analysis_App.py，保存
在当前路径下，点击shift+鼠标右键，在此处点击powershell窗口
输入指令运行：pyinstaller --noconfirm --onefile --windowed --add-data "西南交通大学.png;." --add-data "超导磁浮.jpg;." HTS_Analysis_App.py
在dist文件夹会出现可执行文件
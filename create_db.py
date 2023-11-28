import mysql.connector
from mysql.connector import Error

try:
    connection = mysql.connector.connect(host='localhost',
                                         user='root',
                                         password='')
    user = ''
    password = ''
    
    # Check connection to local mysql server as a root
    if connection.is_connected():
        cursor = connection.cursor()
        # Check whether the database exists, if cancer_uniprotdb does not exist, create one
        cursor.execute("CREATE DATABASE IF NOT EXISTS cancer_uniprotdb")
        cursor.execute("SHOW DATABASES")
        record = cursor.fetchall()
        print("Available databases: ", record)
        # cursor.close()
        # connection.close()

    try:
        # reconnect to local mysql server as a user
        connection = mysql.connector.connect(host='localhost',
                                            user=user,
                                            password=password,
                                            database='cancer_uniprotdb')
        
    except:
        # Grant access to cancer_uniprotdb to user
        cursor.execute("GRANT ALL PRIVILEGES ON cancer_uniprotdb.* TO '{}'@'localhost'".format(user))
        print("Access to cancer_uniprotdb granted for the user")
        # reconnect to local mysql server as a user again
        connection = mysql.connector.connect(host='localhost',
                                            user=user,
                                            password=password,
                                            database='cancer_uniprotdb')
    
    if connection.is_connected():
        print("Connected to database 'cancer_uniprotdb'")
        cursor = connection.cursor()
        cursor.execute("SHOW DATABASES")
        record = cursor.fetchall()

except Error as e:
    print("Error while connecting to MySQL", e)

# finally:
#     print(connection.is_connected())
#     if connection.is_connected():
#         cursor.close()
#         connection.close()
#         print("MySQL connection is closed")